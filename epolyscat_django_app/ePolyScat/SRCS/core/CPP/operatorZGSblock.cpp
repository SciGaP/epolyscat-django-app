// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorZGSblock.h"

#include "tools.h"
#include "useMatrix.h" // this includes the Eigen, which the compiler otherwise does not find

#include "operatorZG.h"

using namespace std;
#include "eigenNames.h"


OperatorZGSblock::OperatorZGSblock(const UseMatrix * Mat, string Kind)
    :OperatorFloor(Mat->rows(),Mat->cols(),"ZGSblock")
{

    // non-zero subblock
    for(unsigned int i=0;i<Mat->rows();i++){if(not Mat->row(i).isZero())rowSub.push_back(i);}
    for(unsigned int i=0;i<Mat->cols();i++){if(not Mat->col(i).isZero())colSub.push_back(i);}

    string hash=Kind+"sparse";
    vector<complex<double> > mDat(Mat->extractSubmatrix(rowSub,colSub));
    sub=new OperatorZG(mDat,rowSub.size(),colSub.size(),hash);

    oNorm=0;
    for(unsigned int k=0;k<mDat.size();k++)oNorm=max(oNorm,abs(mDat[k]));
}
void OperatorZGSblock::pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const{
    Buf.push_back(complex<double>(_rows,_cols));
    Buf.push_back(complex<double>(rowSub.size(),colSub.size()));
    for(unsigned int k=0;k<std::min(rowSub.size(),colSub.size());k++)Buf.push_back(complex<double>(rowSub[k],colSub[k]));
    if(rowSub.size() > colSub.size()) {
        for(unsigned int k = colSub.size(); k < rowSub.size(); ++k)
            Buf.push_back(complex<double>(rowSub[k],0.));
    }
    else if(rowSub.size() < colSub.size()) {
        for(unsigned int k = rowSub.size(); k < colSub.size(); ++k)
            Buf.push_back(complex<double>(0.,colSub[k]));
    }
    Buf.insert(Buf.end(),sub->dat->begin(),sub->dat->end());
    packBasic(Info, Buf);
}
OperatorZGSblock::OperatorZGSblock(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf)
    :OperatorFloor((unsigned int) Buf[0].real(), (unsigned int) Buf[0].imag(), "ZGSblock"){
    auto rowSubSize = (unsigned int) Buf[1].real();
    auto colSubSize = (unsigned int) Buf[1].imag();
    for(unsigned int k=0;k<(unsigned int) rowSubSize;k++)rowSub.push_back((unsigned int)Buf[2+k].real());
    for(unsigned int k = 0; k < colSubSize; ++k)colSub.push_back((unsigned int)Buf[2+k].imag());

    oNorm=DBL_MAX;// dummy norm
    sub=new OperatorZG(vector<complex<double> >(Buf.begin()+std::max(rowSub.size(),colSub.size())+2,Buf.begin()+Info[3]),
            rowSub.size(),colSub.size(),"ZGSblock:"+tools::str(Info[1])+"x"+tools::str(Info[2]));
}

OperatorZGSblock::OperatorZGSblock(const OperatorZGSblock& other)
    : OperatorFloor(other), sub(new OperatorZG(*other.sub)), colSub(other.colSub), rowSub(other.rowSub)
{ }

void OperatorZGSblock::axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
          const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const
{
    vector<complex<double> >y(rowSub.size()),x(colSub.size());
    for(unsigned int k=0;k<colSub.size();k++)x[k]=X[colSub[k]];
    for(unsigned int k=0;k<rowSub.size();k++)y[k]=Y[rowSub[k]];
    sub->axpy(Alfa,x.data(),x.size(),Beta,y.data(),y.size());
    for(unsigned int k=0;k<rowSub.size();k++)Y[rowSub[k]]=y[k];
    if(SizX and SizY)return;
}

void OperatorZGSblock::axpyTranspose(std::complex<double> Alfa, const std::vector<std::complex<double> > &X, std::complex<double> Beta, std::vector<std::complex<double> > &Y) const
{
    vector<complex<double> >x(rowSub.size()),y(colSub.size());
    for(unsigned int k=0;k<rowSub.size();k++)x[k]=X[rowSub[k]];
    for(unsigned int k=0;k<colSub.size();k++)y[k]=Y[colSub[k]];
    sub->axpyTranspose(Alfa,x,Beta,y);
    for(unsigned int k=0;k<colSub.size();k++)Y[colSub[k]]=y[k];
}


