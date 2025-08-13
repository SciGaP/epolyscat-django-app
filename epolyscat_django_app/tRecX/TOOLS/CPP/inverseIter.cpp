// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "inverseIter.h"
#include "useMatrix.h"
#include "abort.h"
#include "tools.h"
#include "printOutput.h"

using namespace std;

void InverseIter::eigen(UseMatrix &Eval, UseMatrix &Evec){

    complex<double> eguess=Eval(0);

    if(Evec.size()!=0){
        if((Evec.rows()!=mat.cols() or Evec.cols()!=Eval.size()))
            ABORT("guess vectors do not match number of Eval or dimension of Mat");
        rvec.resize(Evec.size());
        for(unsigned int k=0;k<rvec.size();k++)rvec[k]=*(Evec.data()+k);
    }

    aMinusEgM=ovr;
    aMinusEgM*=-eguess;
    aMinusEgM+=mat;
    Arpack::eigenIter(Eval.size(),"SmallReal",Evec.size()!=0,true);

    for(unsigned int k=0;k<eval.size();k++)cout<<"eval "<<eval[k]<<"..."<<1./eval[k]<<endl;

    Evec=UseMatrix::UseMap(rvec.data(),mat.cols(),Eval.size());

   if(ovr.size()!=0)mat.eigenOrthonormalize(Eval,Evec,ovr);
   else             mat.eigenOrthonormalize(Eval,Evec,UseMatrix::Identity(mat.rows(),mat.cols()));

    (Evec.transpose()*Evec).print("one?");

    for(unsigned int k=0;k<Evec.cols();k++)
        Eval(k)=(Evec.col(k).transpose()*mat*Evec.col(k))(0)/
                (Evec.col(k).transpose()*ovr*Evec.col(k))(0);

}

void InverseIter::apply(const std::complex<double> *X, std::complex<double> *Y){

    for(unsigned int k=0;k<xTemp.size();k++)*(xTemp.cdata+k)=*(X+k);
    yTemp=ovr*xTemp;
    aMinusEgM.solve(yTemp);
    for(unsigned int k=0;k<xTemp.size();k++)*(Y+k)=*(xTemp.cdata+k);
}

void InverseIter::test(){
    PrintOutput::title("InverseIteration test");
    // size of test
    unsigned int ndim=20;

    // complex symmetric matrix
    UseMatrix mat=UseMatrix::Random(ndim,ndim);
    mat+=UseMatrix::Random(ndim,ndim)*complex<double>(0.,1.);
    UseMatrix ovr=UseMatrix::Random(ndim,ndim);
    ovr=ovr.transpose()*ovr;

    // standard eigensolver
    UseMatrix eval0,evec0;
    mat.eigen(eval0,evec0,ovr);

    // inverse iteration
    UseMatrix eval,evec;
    eval=UseMatrix::Constant(1,1,0.);
    evec=evec0.block(0,0,evec0.rows(),eval.size());
    InverseIter invit(mat,ovr);
    invit.eigen(eval,evec);

    // locate inverse itertion results in excat
    for (unsigned int k=0;k<eval.size();k++){
        bool found=false;
        for (unsigned int l=0;l<eval0.size();l++){
            found=abs(eval0(k)-eval(k))<1.e-5;
            if(found)break;
        }
        if(not found){
            eval0.transpose().print("all eigenvalues");
            ABORT("could not locate eigenvalue "+tools::str(eval(k).complex()));
        }
        else PrintOutput::message("inverse iteration verified eigenvalue "+tools::str(eval(k).complex()));
    }

}
