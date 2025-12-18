// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorRG.h"

#include "tools.h"
#include "readInput.h"
#include "useMatrix.h" // this includes the Eigen, which the compiler otherwise does not find

using namespace std;
#include "eigenNames.h"


OperatorRG::OperatorRG(const std::vector<std::complex<double> > &Phas, const std::vector<double> &mDat, unsigned int Rows, unsigned int Cols, string Kind)
    :OperatorFloor(Rows,Cols,"RG")
{

    string hash=Kind+tools::str(Rows)+tools::str(Cols);

    // store matrix and transposed
    dat=addReal(hash,mDat);
    vector<double> tp(mDat.size());
    for(unsigned int k=0;k<Rows;k++)
        for(unsigned int l=0;l<Cols;l++)
            tp[l+k*Cols]=mDat[k+l*Rows];
    transDat=addReal(hash+"^T",tp);

    x.resize(2*Cols);
    y.resize(2*Rows);

    // remove Phases=1. and =0.
    idx.clear();
    vector<complex<double> > p;
    for(unsigned int k=0;k<Phas.size();k++){
        if(abs(Phas[k]-1.)>1.e-12 and std::norm(Phas[k])>1.e-20){
            p.push_back(Phas[k]);
            idx.push_back(k);
        }
    }
    if(p.size()==0)phas=0;
    else           phas=addComplex(hash+"Phas",p);

    tempComplex.resize(Cols);
    ABORT("what is going on below this line?");
    _rows=Cols;
    _cols=Rows;
    oNorm=0;
    for(unsigned int k=0;k<mDat.size();k++)oNorm=max(oNorm,abs(mDat[k]));
}
void OperatorRG::pack(vector<int> & Info,vector<std::complex<double> > &Buf) const {
    // pack real data into complex storage
    Buf.clear();
    for(int k=0;k<dat->size()-1;k+=2)Buf.push_back(complex<double>(dat->data()[k],dat->data()[k+1]));
    if(dat->size()%2==1)Buf.push_back(dat->back());
    packBasic(Info,Buf);
}

//void OperatorRG::axpy(std::complex<double> Alfa, const std::vector<std::complex<double> > &X,
//                      std::complex<double> Beta, std::vector<std::complex<double> > &Y) const {
//    axpy(Alfa,X.data(),X.size(),Beta,Y.data(),Y.size());
//}

void OperatorRG::axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
                      const std::complex<double> &Beta, std::complex<double>*Y, unsigned int SizY) const {

    if(Alfa==0.){
        scale(Beta,Y,SizY);
        return;
    }

    // multiply by real matrix
    if(tempComplex.size()<SizY)tempComplex.resize(SizY);
    Map<VectorXcd>(tempComplex.data(),SizY).noalias()=Map<MatrixXd>(dat->data(),SizY,SizX)*Map<VectorXcd>(const_cast<complex<double>*>(X),SizX);

    // multiply by phase where needed
    for(unsigned int k=0;k<idx.size();k++)tempComplex.data()[idx[k]]*=phas->data()[k];

    // add into final vector
    if(Beta!=0.){
        scale(Beta,Y,SizY);
        if(Alfa==1.)
            for(complex<double> *y=Y,*a=tempComplex.data();y!=Y+SizY;y++,a++)*y+=*a;
        else
            for(complex<double> *y=Y,*a=tempComplex.data();y!=Y+SizY;y++,a++)*y+=*a*Alfa;
    }
    else {
        if(Alfa==1.)
            for(complex<double> *y=Y,*a=tempComplex.data();y!=Y+SizY;y++,a++)*y=*a;
        else
            for(complex<double> *y=Y,*a=tempComplex.data();y!=Y+SizY;y++,a++)*y=*a*Alfa;
    }

}

void OperatorRG::axpyTranspose(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
                      const std::complex<double> &Beta, std::complex<double>*Y, unsigned int SizY) const {

    if(Alfa==0.){
        scale(Beta,Y,SizY);
        return;
    }

    // multiply by real matrix
    if(tempComplex.size()<SizY)tempComplex.resize(SizY);
    Map<VectorXcd>(tempComplex.data(),SizY).noalias()=Map<VectorXcd>(const_cast<complex<double>*>(X),SizX)*Map<MatrixXd>(dat->data(),SizY,SizX);

    // multiply by phase where needed
    for(unsigned int k=0;k<idx.size();k++)tempComplex.data()[idx[k]]*=phas->data()[k];

    // add into final vector
    if(Beta!=0.){
        scale(Beta,Y,SizY);
        if(Alfa==1.)
            for(complex<double> *y=Y,*a=tempComplex.data();y!=Y+SizY;y++,a++)*y+=*a;
        else
            for(complex<double> *y=Y,*a=tempComplex.data();y!=Y+SizY;y++,a++)*y+=*a*Alfa;
    }
    else {
        if(Alfa==1.)
            for(complex<double> *y=Y,*a=tempComplex.data();y!=Y+SizY;y++,a++)*y=*a;
        else
            for(complex<double> *y=Y,*a=tempComplex.data();y!=Y+SizY;y++,a++)*y=*a*Alfa;
    }

}
