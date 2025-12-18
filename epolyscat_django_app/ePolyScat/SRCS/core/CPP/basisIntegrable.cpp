// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisIntegrable.h"

#include "readInput.h"
#include "basisAbstract.h"
#include "useMatrix.h"
#include "basisMat1D.h"
#include "printOutput.h"
#include "tools.h"

double BasisIntegrable::infty=DBL_MAX;


std::string BasisIntegrable::strDefinition() const{
    return Str(BasisAbstract::strDefinition(),",")
            +order()+lowBound()+upBound()+("["+_comSca->strDefinition()+"]")+_jac->kind();
}

std::string BasisIntegrable::str(int Level) const{
    std::string res=BasisAbstract::str(Level);
    return Str(name(),"")+" "+size()+", "+order()+" ["+_lowBound+","+_upBound+"]";
}

BasisIntegrable::BasisIntegrable(std::string Def){
    std::vector<std::string> def=tools::splitString(Def,',',"[(","])");
    _name=tools::cropString(def[0]);
    _lowBound=tools::string_to_double(tools::cropString(def[3]));
    _upBound=tools::string_to_double(tools::cropString(def[4]));
    _comSca.reset(new ComplexScaling(tools::stringInBetween(def[5],"[","]")));
    _jac.reset(Jacobian::factory(def[6],_comSca->etaX(_lowBound/2.+_upBound/2.)));
}

void BasisIntegrable::valDer(const UseMatrix & X,UseMatrix & Val, UseMatrix & Der, bool ZeroOutside) const{
    std::vector<std::complex<double> >v,d,x;
    for(size_t k=0;k<X.size();k++)x.push_back(X.data()[k]);
    valDer(x,v,d,ZeroOutside);
    Val=UseMatrix::UseMap(v.data(),x.size(),v.size()/x.size());
    Der=UseMatrix::UseMap(d.data(),x.size(),d.size()/x.size());
}
void BasisIntegrable::valDerD(const std::vector<double > & X,
                              std::vector<std::complex<double> > & Val,
                              std::vector<std::complex<double> > & Der, bool ZeroOutside) const{
    std::vector<std::complex<double> >z;
    for(double d: X)z.push_back(d);
    valDer(z,Val,Der,ZeroOutside);
}

UseMatrix BasisIntegrable::val(const std::vector<double> &Grid) const{
    UseMatrix grid(Grid.size(),1);
    for(size_t k=0;k<grid.size();k++)grid(k)=Grid[k];
    return val(grid);
}
UseMatrix BasisIntegrable::val(const UseMatrix& Grid, bool ZeroOutside) const{
    UseMatrix v,d;
    valDer(Grid,v,d,ZeroOutside);
    return v;
}

// this is really bad - we should go directly through the std::vector
std::vector<std::complex<double> > BasisIntegrable::val(double X) const {
    UseMatrix v,d,x(1,1);
    x(0,0)=X;
    valDer(x,v,d);
    return std::vector<std::complex<double> >(v.data(),v.data()+v.size());
}

// this is really bad - we should go directly through the std::vector
std::vector<std::complex<double> > BasisIntegrable::der(double X) const {
    UseMatrix v,d,x(1,1);
    x(0,0)=X;
    valDer(x,v,d);
    return std::vector<std::complex<double> >(d.data(),d.data()+d.size());
}

UseMatrix BasisIntegrable::der(const UseMatrix &Coordinates, bool ZeroOutside) const{
    UseMatrix v,d;
    valDer(Coordinates,v,d,ZeroOutside);
    return d;
}

std::complex<double> BasisIntegrable::eta() const {
    return _comSca->etaX(0.5*(upBound()+lowBound()));
}

void BasisIntegrable::quadRule(int N, UseMatrix &QuadX, UseMatrix &QuadW) const {
    std::vector<double> x,w;
    quadRule(N,x,w);
    if(x.size()!=size_t(N))DEVABORT(Sstr+"quadrature size"+x.size()+"!= N="+N+"for basis"+str());
    QuadX=UseMatrix(N,1);
    QuadW=UseMatrix(N,1);
    for(size_t k=0;k<x.size();k++){
        QuadX(k)=x[k];
        QuadW(k)=w[k];
    }
}

//void BasisIntegrable::quadRule(int N, std::vector<double> &QuadX, std::vector<double> &QuadW) const{
//    UseMatrix x,w;
//    quadRule(N,x,w);
//    for(size_t i=0;i<x.size();i++){
//        QuadX.push_back(x(i).real());
//        QuadW.push_back(w(i).real());
//    }
//}

std::vector<std::vector<std::complex<double> > > BasisIntegrable::transZeroValDer(double Q) const{
    // evaluate values and derivatives at Q
    std::vector<std::complex<double> >vv,dd;
    valDer(std::vector<std::complex<double> >(1,Q),vv,dd,false);

    // single non-zero function
    size_t iNzVal(0);
    double nzVal=0.;
    for(size_t k=0;k<vv.size();k++){
        if(vv[k].imag()!=0.)DEVABORT("not for complex basis: "+str());
        if(abs(vv[k])>1.e-12){
            if(nzVal!=0.)DEVABORT(Str("only evaluate at node of DVR basis, values at Q are:")+vv);
            nzVal=std::abs(vv[k]);
            iNzVal=k;
        }
    }

    // largest in magnitude derivative
    size_t iNzDer(0);
    double nzDer=0.;
    for(size_t k=0;k<dd.size();k++){
        if(k!=iNzVal and nzDer<abs(dd[k]))iNzDer=k;
        nzDer=std::max(nzDer,std::abs(dd[k]));
    }
    if(nzDer==0.)DEVABORT("only zero derivatives?");

    // initial transformation (see notes Stiffness.pdf)
    UseMatrix dMat=UseMatrix::Identity(size(),size());
    for(size_t col=0;col<size();col++)
        if(col!=iNzDer)dMat(iNzDer,col)=dd[col]/dd[iNzDer];

    // label as last and next-to-last functions
    for(size_t k=0;k<size();k++){
        std::swap(dMat.data()[size()*iNzVal+k],dMat.data()[size()*(size()-1)+k]);
        std::swap(dMat.data()[size()*iNzDer+k],dMat.data()[size()*(size()-2)+k]);
    }

    // get metric and (bottom up) Gram-Schmidt orthonormalize
    BasisMat1D ovr("1",this,this);
    UseMatrix uOvr;
    uOvr=ovr.useMat();
    dMat.gramSchmidt(uOvr,std::vector<unsigned int>(),true);

    // put non-zeros into 0,1
    for(size_t k=0;k<size();k++){
        std::swap(dMat.data()[k],dMat.data()[size()*(size()-1)+k]);
        std::swap(dMat.data()[iNzDer+k],dMat.data()[size()*(size()-2)+k]);
    }

    std::vector<std::vector<std::complex<double > > >wMat(size());
    for(size_t k=0;k<size();k++)
        for(size_t l=0;l<size();l++)
            wMat[k].push_back(dMat(l,k).complex());
    return wMat;
}
void BasisIntegrable::plot(std::string Dir) const{
//    if(dynamic_cast<const BasisSet*>(this))return; // BasisSet is NOT plotted here
    UseMatrix q(201,1),val,der;

    // determine reasonable plot range
    double infR=2.;
    double q0=lowBound()<-DBL_MAX/2?upBound()-infR:lowBound();
    double q1= upBound()> DBL_MAX/2?lowBound()+infR:upBound();
    if(q0> DBL_MAX/2)q0= infR;
    if(q1<-DBL_MAX/2)q0=-infR;

    for(unsigned int k=0;k<q.rows();k++)q(k)=q0+k*(q1-q0)/(q.rows()-1);
    valDer(q,val,der);
    std::ofstream out;
    std::string file=Dir+"Basis:"+name()+"_"+tools::str(lowBound(),2);
    for(auto s: {"[","]",",",":"})file=tools::substringReplaceAll(file,s,"_");
    out.open(file.c_str());
    if (not out.is_open())ABORT("could not open file "+file);
    for(unsigned int k=0;k<val.rows();k++){
        out<<q(k).real();
        for (unsigned int l=0;l<val.cols();l++)
            out<<", "<<val(k,l).real()<<", "<<der(k,l).real();
        out<<std::endl;
    }
}
