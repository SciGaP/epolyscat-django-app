// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coefficientsFunction.h"

#include "index.h"
#include "pulse.h"
#include "coefficients.h"
#include "coefficientsFloor.h"
#include "algebra.h"
#include "tools.h"
#include "recursiveIntegrator.h"
#include "parameters.h"
#include "basisGrid.h"
#include "qtEigenDense.h"



using namespace std;
#include "eigenNames.h"

CoefficientsFunction::Data::Data(Tree<Data> *Up, Index *I){
    if(I->parent()!=0 and I->parent()->hasFloor())return;
    for(unsigned int k=0;k<I->childSize() and not I->hasFloor();k++)childAdd(new Data(this,I->child(k)));
}

void Identity::multiply(const Coefficients *X, Coefficients *Y, double Time){*Y=*X;}

void GridValues::multiply(const Coefficients *X, Coefficients *Y, double Time){

    if(X->isLeaf()){
        if(Y->idx()!=X->idx())ABORT("in- and out-index must be identical");
        // apply on floor level
        if(not X->idx()->basis()->grid())ABORT("only for grid discretizatios, found: "+X->idx()->str());
        if(X->idx()->height()!=1)ABORT("implemented only for 1d, is: "+X->idx()->str());
        for(int k=0;k<X->size();k++)Y->floorData()[k]=const_cast<Coefficients*>(X)->floorData()[k]*X->idx()->basis()->grid()->mesh()[k];
    }
    // descend to floor
    for(int k=0;k<X->childSize();k++)multiply(X->child(k),Y->child(k));
}
VolkovPhase::VolkovPhase(const Index* Idx, string kLevel){
    //  if(D->name!="Phi.Eta.Rn_surf_grid")ABORT("gauge transformation not defined for "+D->name);
    const Index * copy;
    copy = Idx;
    kLev = kLevel;
    PhaseAccumulated = new Coefficients(const_cast<Index *>(copy),1.0);
    data=new Data(0,const_cast<Index *>(copy));
    int no_grid=0;

    for(const Index* s = copy;; s=s->child(0)){
        if(s->basis()->isGrid() and s->axisName().substr(0,4)!="spec"
                and s->axisName().substr(0,1)!="k")
            no_grid++;
        if(s->isBottom())break;
    }
    if(no_grid!=0 and no_grid!=2) ABORT(Str("Too many grid axes, cannot identify"," ")
                                        +tools::str(no_grid)+"\n"+copy->hierarchy()+copy->size());

    if(no_grid==2){
        vector<complex<double> >z(no_grid);
        setup(data,copy,z);
    }

    t0 = Parameters::currentTime();  // IMP this should have been set correctly
}

void VolkovPhase::integrateVectorPotentials(double t){
    static tools::RecursiveIntegrator<std::complex<double> > iA_x(Pulse::iAx, 1.e-10),iA_y(Pulse::iAy, 1.e-10),iA_z(Pulse::iAz, 1.e-10);
    vector<double> boundariesOfTimeIntegration(2);
    boundariesOfTimeIntegration[0]=t0; boundariesOfTimeIntegration[1]=t;
    Int_iAx=iA_x.integrate(boundariesOfTimeIntegration,10);
    Int_iAy=iA_y.integrate(boundariesOfTimeIntegration,10);
    Int_iAz=iA_z.integrate(boundariesOfTimeIntegration,10);

}

void VolkovPhase::setup(Data *Dat, const Index* I, vector<complex<double> > &Z){

    if(Dat->root()==Dat){
        vector<double> kgrid;

        const Index *kIdx=I;
        while(kIdx->axisName().find("k")!=0 and not kIdx->isBottom())
            kIdx=kIdx->child(0);
        if(kIdx->axisName().find("k")!=0)ABORT("did not find k-level in hierarcy: "+I->hierarchy());

        for(size_t k=0;k<kIdx->basis()->grid()->mesh().size();k++)
            kgrid.push_back(kIdx->basis()->grid()->mesh()[k]);
        double tol=(kgrid.back()-kgrid[0])*1.e-12;
        bool eqdist=true;
        for(size_t k=1;k<kgrid.size();k++)
            eqdist=eqdist and (abs(kgrid[k]-kgrid[k-1]-kgrid[0])<tol);
        if(eqdist){
            Dat->data.resize(2);
            Dat->data[0].push_back(kgrid[1]-kgrid[0]);
            Dat->data[0].push_back((kgrid[0]+tol)/(Dat->data[0][0]));
        }
        else if(kgrid.size()>10){
            PrintOutput::warning("SLOW: with option -eGrid k-grid is non-equidistant");
        }
    }

    if(Dat->isLeaf()){
        Dat->data.assign(I->sizeStored(),vector<complex<double> >(3));
        for(unsigned int k=0;k<I->sizeStored();k++){
            Dat->data[k][0]=Z[1];                          // z-component
            Dat->data[k][1]=cos(Z[0])*sqrt(1.-Z[1]*Z[1]);  // x-component
            Dat->data[k][2]=sin(Z[0])*sqrt(1.-Z[1]*Z[1]);  // y-component
        }
    }
    else{
        for(unsigned int k=0;k<I->childSize();k++){
            if(Z.size()==2){
                if(I->basis()->isGrid()){
                    if(I->axisName().find("Phi")!=string::npos){
                        Z[0] = I->basis()->grid()->mesh()[k];
                    }
                    else if(I->axisName().find("Eta")!=string::npos){
                        Z[1] = I->basis()->grid()->mesh()[k];
                    }
                    else if(I->axisName().substr(0,4)=="spec" or I->axisName().substr(0,1)=="k")
                    {} //Do nothing
                    else ABORT("Unknown case "+I->axisName());
                }
            }
            else ABORT("Unknown case");
            setup(Dat->at(k),I->child(k),Z);
        }
    }
}

TIMER(integ,volkov)
TIMER(update,)
TIMER(multip,volkov)
void VolkovPhase::multiply(const Coefficients *X, Coefficients *Y, double Time){
    STARTDEBUG(integ);
    integrateVectorPotentials(Time);  // setup vector potential integrals
    STOPDEBUG(integ);
    STARTDEBUG(update);//here
    if(Time != t0) UpdateVolkovPhase(data,Time,PhaseAccumulated);
    STOPDEBUG(update);
    STARTDEBUG(multip);//here
    multiplyLocal(X,Y,PhaseAccumulated);
    STOPDEBUG(multip);
    t0 = Time;                        // Set t0 to current time
}

void VolkovPhase::consistency(std::string Operator){
    if(Operator.find("-iLaserA"))
        PrintOutput::warning("field-coupling through -iLaserA*[t] may be inconsistent with Volkov:\nH(t) = "
                             +Operator+"\nVolkov phase assumes -Laplace/2 + i A[t].Nabla + ...");
}

void VolkovPhase::UpdateVolkovPhase(Data *Dat, double Time, Coefficients *phA){
    double dTimeHalf=0.5*(Time-t0);
    if(Dat->root()==Dat){
        if(Dat->data.size()>0){
            const Index* s=phA->firstLeaf()->idx();
            if(not s->basis()->grid())s=s->descend();//HACK
            Dat->data[1].clear();
            for(unsigned int i = 0; i < s->basis()->grid()->size(); i++)
                Dat->data[1].push_back(s->basis()->grid()->mesh()[i]);
            for(size_t k=0;k<Dat->data[1].size();k++)
                Dat->data[1][k]=exp(complex<double>(0.,pow(Dat->data[1][k].real(),2)*dTimeHalf));
        }
    }
    if(not phA->idx()->hasFloor())
        for(unsigned int k=0;k<phA->childSize();k++){
            UpdateVolkovPhase(Dat->at(k),Time,phA->child(k));
        }
    else{
        // Number of floor levels
        int floorlevels=0;
        for(const Index* s=phA->idx();;s=s->child(0)){floorlevels++;if(s->isBottom())break;}

        // Volkov Phase
        for(const Index* s=phA->idx();;s=s->child(0)){
            if(s->axisName()==kLev){  // Momentum index, multiply by a volkov phase
                vector<complex<double> > K;
                for(double p: s->basis()->grid()->mesh())K.push_back(p);
                complex<double > kx = Int_iAx;
                complex<double > ky = Int_iAy;
                complex<double > kz = Int_iAz;
                complex<double> dExp(0.),expVolk;
                if(Dat->data.size()!=0){ // Spherical coordinates
                    // note: data only contains unit vectors, but floor is at fixed phi,eta
                    kz *= Dat->data[0][0];
                    kx *= Dat->data[0][1];
                    ky *= Dat->data[0][2];
                    if(Dat->root()->data.size()>0){
                        dExp=exp(-(kx+ky+kz)*Dat->root()->data[0][0]);
                        expVolk=pow(dExp,int(Dat->root()->data[0][1].real()*1.00000001));
                    }
                }
                if(floorlevels==1){
                    if(K.size() != phA->size()) DEVABORT("K size doesn't match floor size");
                    if(dExp==0.){
                        for(unsigned int i=0;i<phA->size();i++){
                            // NOTE: exponentiation can be avoided if equistant K-grid is used
                            phA->data()[i] *=
                                    exp((complex<double>(0.,K[i].real()*dTimeHalf)-kx-ky-kz)*K[i].real() );
                        }
                    }
                    else {
                        complex<double>*expKsq=Dat->root()->data[1].data();
                        complex<double>*phAdat=phA->data();
                        complex<double>*phAend=phAdat+phA->size();
                        for(;phAdat!=phAend;phAdat++,expKsq++){
                            *phAdat*=expVolk*(*expKsq);
                            expVolk*=dExp;
                        }
                    }
                }
                else if(floorlevels==2){

                    int otherfloorsize= phA->size()/K.size();

                    if(s->depthInFloor()==0){  // k is the first level
                        for(unsigned int i=0;i<K.size();i++){

                            complex<double > K2 = complex<double>(0., 0.5)*pow(K[i].real(), 2)*(Time-t0);
                            complex<double > kx = Int_iAx;
                            complex<double > ky = Int_iAy;
                            complex<double > kz = Int_iAz;

                            if(Dat->data.size()!=0){ // Spherical coordinates
                                kz *= Dat->data[i][0]; kx *= Dat->data[i][1]; ky *= Dat->data[i][2];
                            }

                            Map<VectorXcd>(phA->floorData()+i*otherfloorsize,otherfloorsize) *=
                                    exp(K2 - K[i].real()*(kx+ky+kz));
                        }
                    }
                    else if(s->depthInFloor()==1){  // go through strides
                        for(unsigned int i=0;i<K.size();i++){

                            complex<double > K2 = complex<double>(0., 0.5)*pow(K[i].real(), 2)*(Time-t0);
                            complex<double > kx = Int_iAx;
                            complex<double > ky = Int_iAy;
                            complex<double > kz = Int_iAz;

                            if(Dat->data.size()!=0){ // Spherical coordinates
                                kz *= Dat->data[i][0]; kx *= Dat->data[i][1]; ky *= Dat->data[i][2];
                            }
                            Map<VectorXcd,0,Eigen::InnerStride<> >(phA->floorData()+i,otherfloorsize,Eigen::InnerStride<>(K.size())) *=
                                    exp(K2 - K[i].real()*(kx+ky+kz));
                        }

                    }
                    else DEVABORT("Smothing Wrong");
                }
                else DEVABORT("Case not implemented");
            }
            if(s->isBottom())break;
        }
    }
}

void VolkovPhase::multiplyLocal(const Coefficients *X, Coefficients *Y, Coefficients *phA){
    if(not phA->idx()->hasFloor())
        for(unsigned int k=0;k<X->childSize();k++)
            multiplyLocal(X->child(k),Y->child(k),phA->child(k));
    else
        Map<ArrayXcd>(Y->floorData(),Y->size()) =
                Map<ArrayXcd>(phA->floorData(),phA->size())*Map<ArrayXcd>(X->floorData(),X->size());
}
