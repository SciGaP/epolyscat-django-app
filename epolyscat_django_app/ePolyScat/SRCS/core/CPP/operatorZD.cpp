// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorZD.h"

#include "tools.h"
#include "useMatrix.h" // this includes the Eigen, which the compiler otherwise does not find

using namespace std;
#include "eigenNames.h"


OperatorZD::OperatorZD(const UseMatrix *Mat, string Kind)
    :OperatorFloor(Mat->rows(),Mat->cols(),"ZD"){
    construct(Mat,Kind);
}

void OperatorZD::construct(const UseMatrix *Mat, string Kind)
{
    string hash=Kind+tools::str(Mat->rows())+tools::str(Mat->cols());
    vector<complex<double> >mDat(min(Mat->rows(),Mat->cols()));
    for(unsigned int k=0;k<mDat.size();k++)mDat[k]=Mat->operator()(k,k);
    dat=addComplex(hash,mDat);
    _rows=Mat->rows();
    _cols=Mat->cols();
    oNorm=Mat->maxAbsVal();
}

OperatorZD::OperatorZD(const Eigen::MatrixXcd &Mat, string Kind)
    :OperatorFloor(Mat.rows(),Mat.cols(),"ZD"){
    Eigen::VectorXcd vec(Mat.cols());
    for(int k=0;k<vec.size();k++)vec(k)=Mat(k,k);
    construct(vec,Kind);
}
OperatorZD::OperatorZD(const std::complex<double> * Diag, int Size)
    :OperatorFloor(Size,Size,"ZD"){
    construct(Eigen::Map<Eigen::VectorXcd>(const_cast<std::complex<double>*>(Diag),Size),"OvrDVR");
}

void OperatorZD::construct(const Eigen::VectorXcd & Mat, string Kind)
{
    string hash=Kind+tools::str(Mat.size());
    vector<complex<double> >mDat(Mat.size());
    for(unsigned int k=0;k<mDat.size();k++)mDat[k]=Mat(k);
    dat=addComplex(hash,mDat);
    _rows=Mat.size();
    _cols=Mat.size();
    oNorm=Mat.lpNorm<Eigen::Infinity>();
}

void OperatorZD::pack(vector<int> & Info,vector<std::complex<double> > &Buf) const {
    Buf.insert(Buf.end(),dat->begin(),dat->end());
    packBasic(Info,Buf);
}
void OperatorZD::axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
                      const std::complex<double> &Beta,       std::complex<double> *Y, unsigned int SizY) const {

    if(Alfa==0.){
        scale(Beta,Y,SizY);
    } else if(Alfa==1.){
        if(Beta==1.)     for(unsigned int k=0;k<dat->size();k++)Y[k]+=*(dat->data()+k)*X[k];
        else if(Beta==0.)for(unsigned int k=0;k<dat->size();k++)Y[k] =*(dat->data()+k)*X[k];
        else             for(unsigned int k=0;k<dat->size();k++)Y[k]=Beta*Y[k]+(*(dat->data()+k))*X[k];
    } else {
        if(Beta==1.)     for(unsigned int k=0;k<dat->size();k++)Y[k]+=*(dat->data()+k)*X[k]           *Alfa;
        else if(Beta==0.)for(unsigned int k=0;k<dat->size();k++)Y[k]  =*(dat->data()+k)*X[k]          *Alfa;
        else             for(unsigned int k=0;k<dat->size();k++)Y[k]=Beta*Y[k]+(*(dat->data()+k))*X[k]*Alfa;
    }
}

long OperatorZD::applyCount() const{
    return _rows;
}
