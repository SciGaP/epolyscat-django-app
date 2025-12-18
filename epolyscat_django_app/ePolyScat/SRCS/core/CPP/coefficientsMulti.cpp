// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coefficientsMulti.h"

#include "operatorAbstract.h"
#include "useMatrix.h"
#include "printOutput.h"

using namespace std;

CoefficientsMulti::CoefficientsMulti(const Index * Idx,unsigned int Cols,complex<double> Val){
    if(Cols==0)return;
    c.push_back(Coefficients(Idx,Val));
    for(int k=1;k<Cols;k++)c.push_back(c.back());
}

CoefficientsMulti & CoefficientsMulti::apply(std::complex<double> A,const OperatorAbstract & Op, CoefficientsMulti & C){
    if(C.vecs()!=vecs())ABORT("multiplicities do not match");
    for(int k=0;k<c.size();k++)Op.apply(A,C.vec(k),0.,c[k]);
    return *this;
}

UseMatrix CoefficientsMulti::innerProduct(const CoefficientsMulti & Rhs) const{
    UseMatrix p(vecs(),Rhs.vecs());
    for(int j=0;j<Rhs.vecs();j++)
        for(int i=0;i<vecs();i++)
            p(i,j)=c[i].innerProduct(&Rhs.c[j],true);
    return p;
}
CoefficientsMulti & CoefficientsMulti::operator*=(const UseMatrix & Mat){
    CoefficientsMulti res(c[0].idx(),Mat.rows(),0.);
    for(int j=0;j<Mat.cols();j++)
        for(int i=0;i<Mat.rows();i++)
            res.c[j]+=Mat(i,j).complex()*c[i];
    std::swap(res.c,c);
    return *this;
}

CoefficientsMulti & CoefficientsMulti::orthonormalize(const OperatorAbstract &Overlap, bool Pseudo, double EpsRemove) {
    if(c.size()==0)return *this;
    Coefficients *sCk=Overlap.tempLHS();
    int k=0;
    while(k<c.size()){
        Overlap.apply(1.,c[k],0.,*sCk);
        for(int l=0;l<k;l++)c[k].axpy(-c[l].innerProduct(sCk,Pseudo),&c[l]);
        complex<double> a=c[k].innerProduct(sCk,Pseudo);
        if(abs(a)<EpsRemove){
            c.erase(c.begin()+k);
        }
        else {
            if(abs(a)<1.e-12){
                if(abs(a)==0.)ABORT("vectors linearly dependent");
                PrintOutput::warning("near-linearly dependent vectors");
            }
            c[k]*=1./sqrt(a);
            k++;
        }
    }
    return *this;
}
