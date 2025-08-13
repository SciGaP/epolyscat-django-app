// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "eigenSubspace.h"

#include "readInput.h"
#include "printOutput.h"

#include "operatorTree.h"
//#include "operator.h"
#include "discretizationSpectral.h"
#include "useMatrix.h"
#include "operatorDiagonal.h"
#include "algebra.h"

using namespace std;

EigenSubspace::EigenSubspace(ReadInput &Inp){
    double reE,imE;
    int line=1;
    do { // run at least once for input documentation
        Inp.read("EigenSubspace","imGuess",imE,"0.","real part of spectral center",line);
        Inp.read("EigenSubspace","reGuess",reE,Algebra::Infty,"real part of spectral center",line);
        eResolv.push_back(complex<double>(reE,imE));
    }
    while (not Inp.endCategory("EigenSubspace",++line));
    Inp.read("EigenSubspace","tolerance",tolerance,"1.e-7","largest eigenvalue error (relative to radius)");
    eGuess=eResolv;
}

void EigenSubspace::eigen(const std::vector<std::complex<double> > &E0, std::vector<std::complex<double> > & ETarget, vector<Coefficients*> & Evec){
    if(ETarget.size()==0)ETarget=E0;
    if(ETarget.size()!=E0.size())ABORT("if ETarget is specified, size must match E0");
    if(Evec.size()!=0 and Evec.size()!=ETarget.size())ABORT("if Evec is specified, size must match E0");

    setTarget(E0,ETarget,Evec);
    ETarget=iterate();
    for(int k=0;k<c0.vecs();k++){
        if(k<Evec.size())*Evec[k]=c0.vec(k);
        else  Evec.push_back(new Coefficients(c0.vec(k)));
    }
}

void EigenSubspace::print() const {
    PrintOutput::lineItem("EigenSubspace: tolerance",tolerance);
    PrintOutput::paragraph();

}

bool EigenSubspace::converged(const UseMatrix &Smat){
    // number of target energies converged to tolerance

    if(Smat.size()>0){
        UseMatrix sval;
        Smat.eigenValues(sval);
        double amin=DBL_MAX,amax=0.;
        for(int k=0;k<sval.size();k++){
            amin=min(amin,abs(sval(k)));
            amax=max(amax,abs(sval(k)));
        }
        if(amin<amax*1.e-12)return true;
    }

    if(e0[0].size()<3)return false;
    int conN=0;
    for(int k=0;k<e0.size();k++){
        bool conk=true;
        for(int m=e0[k].size()-3;m<e0[k].size();m++)
            for(int l=e0[k].size()-3;l<m;l++){
                if(abs(e0[k][l]-e0[k][m])>tolerance){
                    conk=false;
                }
            }
        if(conk)conN++;
    }
    return conN>=eResolv.size();
}

void EigenSubspace::setup(const Discretization * D, const OperatorAbstract *H, const OperatorAbstract *H0)
{
    // get spectral discretization for H0
    h=H;
    ABORT("disabled at present");
    //specH0=new DiscretizationSpectral(D,H0);
    vector<Coefficients*>evec;
    setTarget(eResolv,eResolv,evec);
}

complex<double> EigenSubspace::eCenter() const{
    complex<double> e=0;
    for(int k=0;k<eResolv.size();k++)e+=eResolv[k];
    e/=int(eResolv.size());
    return e;
}

void EigenSubspace::setTarget(vector<complex<double> > E0, vector<complex<double> > & ETarget, std::vector<Coefficients *> &Evec){
    if(Evec.size()!=0 and Evec.size()!=ETarget.size())ABORT("number of eigenvalues does not match number of guess vectors");
    eGuess=ETarget;
    eResolv=E0;
    if(eResolv.size()==0)eResolv=eGuess;


    // get the sorting of eigenvalues by minimal distance to E0
    vector<complex<double> > eval;
    for(int k=0;k<specH0->spectralOper()->vals();k++)eval.push_back(specH0->spectralOper()->val(k));
    vector<int> srt;
    sortByDistance(E0,UseMatrix(),eval,srt);

    // get starting C0: eigenvectors of H0 within radius around center
    Coefficients * specC=specH0->mapFromParent()->tempLHS();
    specC->treeOrderStorage();
    c0.clear();
    for(int k=0;k<min(eResolv.size(),srt.size());k++){
//        if(Evec.size()!=0){
//            c0.push_back(*Evec[k]);
//        } else {
            specC->setToZero();
            specC->storageData()[srt[k]]=1.;
            c0.push_back(*specH0->mapFromParent()->tempRHS());
            specH0->mapToParent()->apply(1.,*specC,0.,c0.back());
//        }
        specH0->spectralOper()->setFunctionValue(srt[k],0.);
        cout<<"zeroing "<<specH0->spectralOper()->val(srt[k])<<endl;
    }

    // put (H0-eCenter)^-1

    cout<<"center at "<<eCenter()<<" for "<<tools::str(E0,5,", ")<<endl;
    for(int k=eResolv.size();k<srt.size();k++){
        complex<double> res=specH0->spectralOper()->val(srt[k])-eCenter();
        if(abs(res)<1.e-8)res=1.e-8; // avoid accidental over-emphasis
        specH0->spectralOper()->setFunctionValue(srt[k],1./res);
    }
}

CoefficientsMulti & EigenSubspace::applyResolvent(CoefficientsMulti &C) const{
    for(int k=0;k<C.vecs();k++){
        specH0->mapFromParent()->apply(1.,C.vec(k),0.,*specH0->mapFromParent()->tempLHS());
        specH0->spectralOper()->apply(1.,*specH0->mapFromParent()->tempLHS(),0.,*specH0->mapFromParent()->tempLHS());
        specH0->mapToParent()->apply(1.,*specH0->mapFromParent()->tempLHS(),0.,C.vec(k));
    }
    return C;
}


vector<complex<double> > EigenSubspace::iterate(){

    c0.orthonormalize(*h->iIndex->overlap(),true);
    CoefficientsMulti sC(c0);
    CoefficientsMulti hC(c0);
    sC.apply(1.,*h->iIndex->overlap(),c0);
    hC.apply(1.,*h,c0);
    UseMatrix h00=c0.innerProduct(hC);
    e0.assign(c0.vecs(),vector<complex<double> >(1));
    for(int k=0;k<c0.vecs();k++)e0[k][0]=h00(k,k).complex();

    CoefficientsMulti sC0(sC);
    UseMatrix smat,hmat;
    while(not converged(smat)){
        // D = -(H_C0 - S_C0 C0^T H C0)=S_C0 E0 - hC
        CoefficientsMulti c1(sC);
        h00=c0.innerProduct(hC);
        c1*=h00;
        c1-=hC;

        // remove linear dependencies
        c1.orthonormalize(*h->iIndex->overlap(),true);

        // C1 <- (H0-SEg)^-1 C1
        applyResolvent(c1);

        int d0=c0.vecs(),d1=c1.vecs();
        smat=UseMatrix::Constant(d0+d1,d0+d1,0.);
        hmat=UseMatrix::Constant(d0+d1,d0+d1,0.);
        smat.block(0,0,d0,d0)=c0.innerProduct(sC);
        hmat.block(0,0,d0,d0)=h00;

        // H_C1, S_C1
        CoefficientsMulti hC1(h->iIndex,d1),sC1(h->iIndex,d1);
        hC1.apply(1.,*h,c1);
        sC1.apply(1.,*h->iIndex->overlap(),c1);

        // H01 = C1^T H_C0, S01 = C0^T S_C1
        smat.block(0,d0,d0,d1)=c0.innerProduct(sC1);
        hmat.block(0,d0,d0,d1)=c0.innerProduct(hC1);

        smat.block(d0,0,d1,d0)=smat.block(0,d0,d0,d1).transpose();
        hmat.block(d0,0,d1,d0)=hmat.block(0,d0,d0,d1).transpose();

        // H11 = C1^T H_C1, S11 = C1^T S_C1
        smat.block(d0,d0,d1,d1)=c1.innerProduct(sC1);
        hmat.block(d0,d0,d1,d1)=c1.innerProduct(hC1);

        // solve eigenproblem
        UseMatrix subE,subV;
        hmat.eigen(subE,subV,smat);
        if(not hmat.eigenOrthonormalize(subE,subV,smat)){
            (smat-smat.transpose()).print("smat",0);
            (hmat-hmat.transpose()).print("smat",0);
            ABORT("not hermitian or symmetric");
        }

        // select subset sel of eigenvectors by closeness to target energies
        vector<int> idx;
        vector<complex<double> >vals;
        for(int k=0;k<subE.size();k++)vals.push_back(subE(k).complex());

        // get overlaps of all new vectors with target
        UseMatrix ovrTarget=smat.block(0,0,d0,d0+d1);
//        ovrTarget*=subV;
        Eigen::Map<Eigen::MatrixXcd>(ovrTarget.data(),ovrTarget.rows(),ovrTarget.cols())*=Eigen::Map<Eigen::MatrixXcd>(subV.data(),subV.rows(),subV.cols());

        sortByDistance(eGuess,ovrTarget,vals,idx);

        UseMatrix sel0(d0,d0);
        UseMatrix sel1(d1,d0);
        for(int k=0;k<d0;k++){
            sel0.col(k)=subV.block(0, idx[k],d0,1);
            sel1.col(k)=subV.block(d0,idx[k],d1,1);
            e0[k].push_back(subE(idx[k]).complex());
        }

        // compose new c0, H_C0 and S_C0 as linear combinations (H_C0,H_C1) sel  and (S_C0,S_C1) sel
        c0*=sel0;
        c0+=c1*=sel1;

        sC.apply(1.,*h->iIndex->overlap(),c0);
        hC.apply(1.,*h,c0);
    }
    vector<complex<double> > res;
    for(int k=0;k<eResolv.size();k++){
        res.push_back(e0[k].back());
    }
    return res;
}

void EigenSubspace::sortByDistance(const vector<complex<double> > & Centers, const UseMatrix & OTarget, const vector<std::complex<double> > & Values, vector<int> &Sort){
    vector<double> dist;
    Sort.clear();
    for(int k=0;k<Values.size();k++){
        Sort.push_back(k);
//        if(OTarget.size()>0){
//            dist.push_back(abs(OTarget(0,k).complex()-1.));
//            for(int l=1;l<Centers.size();l++)dist.back()=min(dist.back(),abs(OTarget(l,k).complex()-1.));
//        } else {
            dist.push_back(abs(Values[k]-Centers[0]));
            for(int l=1;l<Centers.size();l++)dist.back()=min(dist.back(),abs(Values[k]-Centers[l]));
//        }
    }
    tools::sortByKey(dist,Sort);
}


void EigenSubspace::write(const string &File){

    PrintOutput::title("EIGENVALUES");
    PrintOutput::newRow();
    PrintOutput::rowItem("real");
    PrintOutput::rowItem("imag");
    for (unsigned int k=0;k<eResolv.size();k++){
        PrintOutput::newRow();
        PrintOutput::rowItem(e0[k].back().real(),12);
        PrintOutput::rowItem(e0[k].back().imag(),12);
    }

    ofstream eig;
    eig.open(File.c_str());
    for (unsigned int n=0;n<eResolv.size();n++){
        eig<<setprecision(10)<<real(e0[n].back())<<", "<<imag(e0[n].back())<<endl;
    }        PrintOutput::title("EIGENVALUES");


    PrintOutput::message(tools::str(eResolv.size())+" eigenvalues on "+ReadInput::main.output()+"eig",0,true);

}
