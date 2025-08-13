// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorZtimesId.h"

#include "tools.h"
#include "useMatrix.h" // this includes the Eigen, which the compiler otherwise does not find

using namespace std;
using namespace Eigen;

OperatorZtimesId::OperatorZtimesId(std::complex<double> Coeff, int Rows, string Kind)
    :OperatorFloor(Rows,Rows,"ZtimesId")
{
    string hash=Kind+tools::str(Rows);
    vector<complex<double> >mDat;
    mDat.push_back(Coeff);
    dat=addComplex(hash,mDat);
    _rows=Rows;
    _cols=Rows;
    oNorm=std::abs(Coeff);
}
void OperatorZtimesId::pack(vector<int> & Info,vector<std::complex<double> > &Buf) const {
    Buf.insert(Buf.end(),dat->begin(),dat->end());
    packBasic(Info,Buf);
}
void OperatorZtimesId::axpy(const std::complex<double>& Alfa, const std::complex<double>* X, unsigned int SizX,
                            const std::complex<double>& Beta, std::complex<double>* Y, unsigned int SizY) const{

    if(Alfa==0.){
        scale(Beta,Y,SizY);
    } else if(Alfa==1.){
        if(Beta==1.)     for(unsigned int k=0;k<_rows;k++)Y[k]+=*(dat->data())*X[k];
        else if(Beta==0.)for(unsigned int k=0;k<_rows;k++)Y[k] =*(dat->data())*X[k];
        else             for(unsigned int k=0;k<_rows;k++)Y[k] =Beta*Y[k]+(*(dat->data()))*X[k];
    } else {
        if(Beta==1.)     for(unsigned int k=0;k<_rows;k++)Y[k]+=*(dat->data())*X[k]           *Alfa;
        else if(Beta==0.)for(unsigned int k=0;k<_rows;k++)Y[k] =*(dat->data())*X[k]          *Alfa;
        else             for(unsigned int k=0;k<_rows;k++)Y[k] =Beta*Y[k]+(*(dat->data()))*X[k]*Alfa;
    }
}

long OperatorZtimesId::applyCount() const{
    return _rows;
}
