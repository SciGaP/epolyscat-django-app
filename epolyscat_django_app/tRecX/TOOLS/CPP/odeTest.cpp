// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "odeTest.h"

#include "odeSIA.h"
#include "odeCorSIA.h"
#include "odeRK4.h"
#include "vectorComplex.h"
#include "linSpaceMap.h"
#include "vectorComplex.h"
#include "arnoldi.h"
#include "useMatrix.h"

using namespace std;

TestDer::TestDer(const VectorComplex & Dat, const VectorComplex & Ht, const VectorComplex & Vec):dat(Dat),ht(Ht),vec(Vec),tCurr(0.){}

void TestDer::update(double Time, const std::vector<std::complex<double>> & Parameters){
    dat.axpy(Time-tCurr,ht,1.);
    tCurr=Time;
}
void TestDer::update(double Time, const VectorComplex* Vec){
    dat.axpy(Time-tCurr,ht,1.);
    tCurr=Time;
}

void TestDer::apply(std::complex<double> A, const VectorComplex &Vec, std::complex<double> B, VectorComplex &Y) const{
    Y*=B;
    const std::complex<double> *m=dat.data();
    for (unsigned int l=0;l<Vec.size();l++)
        for (unsigned int k=0;k<Y.size();k++,m++)
            Y[k]+=*m*Vec[l]*A;
}
const VectorComplex & TestDer::lhsVector() const {return vec;}

/// test OdeSIA (short iterative arnoldi)
void OdeTest::SIA(){

    unsigned int dim=500;

    VectorComplex dat(dim*dim),ht(dim*dim),vec(dim);
    for(unsigned int k=0;k<dat.size();k++){
        dat[k]=std::complex<double>((std::rand() % 10000)*1.e-4,(std::rand() % 10000)*1.e-4);
        ht[k]=std::complex<double>((std::rand() % 10000)*1.e-4,(std::rand() % 10000)*1.e-4);
    }

    VectorComplex ovr(dim);
    for(unsigned int k=0;k<ovr.size();k++)ovr[k]=0.1+abs(std::rand() %1000*1e-4);

    //    // make hermitian
    //    for(unsigned int i=0;i<dim;i++){
    //        dat[i+i*dim].imag()=0.;
    //        ht[i+i*dim].imag()=0.;
    //        for(unsigned int j=0;j<i;j++){
    //            dat[i+j*dim]=conj(dat[j+i*dim]);
    //            ht[i+j*dim]=conj(ht[j+i*dim]);
    //        }
    //    }

    // left multiply by "inverse overlap" and -i
    complex<double> eta=exp(std::complex<double>(0.,-0.5));
    dat*=eta;
    ht *=sqrt(eta);
    for(unsigned int i=0,ij=0;i<dim;i++){
        for(unsigned int j=0;j<dim;j++,ij++){
            dat[ij]*=complex<double>(0.,-1.)/ovr[j];
            ht[ij] *=complex<double>(0.,-1.)/ovr[j];
        }
    }

    // random initial vector
    for(unsigned int k=0;k<vec.size();k++)
        vec[k]=std::complex<double>((std::rand() % 10000)*1.e-4,(std::rand() % 10000)*1.e-4);

    TestDer der(dat,ht,VectorComplex(dim));
    Arnoldi<TestDer,VectorComplex> arn(der,vec);
    UseMatrix val0,id=UseMatrix::Identity(dim,dim);
    unsigned int dim0=4;
    for(unsigned int it=0;it<1;it++){
        arn.reset(vec);
        cout<<"extend from "<<arn.krylovDim();
        arn.extend(dim0);
        cout<<" to "<<arn.krylovDim();
        arn.extend(dim);
        cout<<" to "<<arn.krylovDim()<<endl;

        // first compare eigenvalues
        UseMatrix val1;
        UseMatrix::UseMap(dat.data(),dim,dim).eigenValues(val0,id);
        UseMatrix::UseMap(arn.matrix().data(),dim,dim).eigenValues(val1,id);

        // sort by increasing real part
        UseMatrix val00(val0),val11(val1);
        std::sort(val00.data(),val00.data()+val00.rows(),tools::lessReal);
        std::sort(val11.data(),val11.data()+val11.rows(),tools::lessReal);

        if(not (val00-val11).isZero(1.e-9)){
            (val00-val11).transpose().show("diff: "+tools::str((val00-val11).maxAbsVal()));
        }
        else cout<<"Arnoldi matrix OK"<<endl;
        cout<<"extremal eigenvalues "<<val00(0).complex()<<" "<<val00(val00.rows()-1).complex()<<endl;
    }

    // time-propagate exact
    UseMatrix vecR,vecL;
    double time=0.004;
    dat.axpy(time/2.,ht,1.);

    UseMatrix h;
    h=UseMatrix::UseMap(dat.data(),dim,dim);
    h.eigen(val0,vecR,id);
    // check eigenvectors
    UseMatrix valMat=UseMatrix::Zero(dim,dim);
    for(unsigned int k=0.;k<dim;k++)valMat(k,k)=val0(k);
    cout<<"max eigenvector err "<<(vecR*valMat-h*vecR).maxAbsVal()<<endl;

    vecL=vecR.inverse().adjoint();
    cout<<"identity L^dagger R: "<<(vecL.adjoint()*vecR).isIdentity(1.e-12)<<endl;
    UseMatrix vec0,vecX,vecA;
    vec0=UseMatrix::UseMap(vec.data(),dim,1);
    vecX=vecL.adjoint()*vec0;
    for(unsigned int k=0;k<val0.size();k++)vecX(k)*=exp(time*val0(k).complex());
    vecX=vecR*vecX;


    // time-propagate Arnoldi
    OdeSIA<TestDer,VectorComplex> sia(&der,dim,1.e-12);
    //    OdeCorSIA<TestDer,VectorComplex> sia(&der,dim,1.e-12);
    sia.step(vec,0.,time);
    vecA=UseMatrix::UseMap(vec.data(),dim,1);
    cout<<"\n"+sia.name()+" - max diff "<<(vecX-vecA).maxAbsVal()<<" at "<<sia.currentKrylov()<<" (max="<<sia.maxKrylov()<<")"<<endl;

    // step-size dependence of error estimate
    bool first=true;
    double quadErr=0.;
    VectorComplex err(der.lhsVector());

    cout<<"\n"+sia.name()+" consistency: "+tools::str(sia.consistencyOrder())
          +"\n Step,  Estimate, "+tools::str(sia.consistencyOrder()+1)+"th power growth"<<endl;
    double tcur=time/64;
    first=true;
    VectorComplex vSia;
    while(tcur<time*2.){
        vSia=vec; // always use same vector
        sia.stepError(vSia,0,tcur,err);
        if(first)quadErr=err.norm();
        else quadErr*=std::pow(2,sia.consistencyOrder()+1);
        cout<<tcur<<", "<<err.norm()<<", "<<quadErr<<" "<<sia.nApplyStep()<<endl;
        tcur*=2.;
        first=false;
    }


    OdeRK4<TestDer,VectorComplex> rk4(&der);
    cout<<"\n"+rk4.name()+" consistency: "+tools::str(rk4.consistencyOrder())
          +"\n Step,  Estimate, "+tools::str(rk4.consistencyOrder()+1)+"th power growth"<<endl;
    tcur=time/256;
    first=true;
    VectorComplex vRk4(vec);
    while(tcur<time*2.){
        vRk4=vec;
        vSia=vRk4;
        sia.stepError(vSia,0,tcur,err);
        rk4.stepError(vRk4,0,tcur,err);
        if(first)quadErr=err.norm();
        else quadErr*=std::pow(2,rk4.consistencyOrder()+1);
        cout<<tcur<<", "<<err.norm()<<", "<<quadErr<<" "<<rk4.nApplyStep()<<" "<<(vSia-=vRk4).norm()<<endl;
        tcur*=2.;
        first=false;
    }

    // compare
    vRk4=vec;
    vSia=vRk4;
    time=0.002;
    unsigned int n=20;
    for(unsigned int k=0;k<n;k++){
        rk4.stepError(vRk4,k*time/n,time/n,err);
    }
    sia.stepError(vSia,0,time,err);
    cout<<"\n multistep RK4 vs. sia: "<<time<<", Err(sia)="<<err.norm()<<" Sia-RK4="<<(vSia-=vRk4).norm()<<" apply(Sia)="<<sia.nApplyStep()<<endl;






}
