// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorTensorProduct.h"

#include "tools.h"
#include "useMatrix.h" // this includes the Eigen, which the compiler otherwise does not find

using namespace std;
#include "eigenNames.h"


OperatorTensorProduct::OperatorTensorProduct(const std::vector<int> &Info, const std::vector<std::complex<double> >&Buf,std::string Kind)
    :OperatorFloor(Kind){
    unpackBasic(Info,Buf);
    subRows.push_back((unsigned int)Buf[0].real());
    subCols.push_back((unsigned int)Buf[0].imag());
    subRows.push_back(_rows/subRows[0]);
    subCols.push_back(_cols/subCols[0]);
}

OperatorTensorProduct::OperatorTensorProduct(std::vector<const UseMatrix*> Dat, string Kind)
    :OperatorFloor(0,0,Kind)
{
    oNorm=1.;
    _rows=1;
    _cols=1;
    for(unsigned int i=0;i<Dat.size();i++){
        subRows.push_back(Dat[i]->rows());
        subCols.push_back(Dat[i]->cols());
        _rows*=subRows.back();
        _cols*=subCols.back();
        oNorm*=Dat[i]->maxAbsVal();
    }
}
void OperatorTensorProduct::pack(std::vector<int> &Info, std::vector<complex<double> >&Buf) const{
    Buf.push_back(complex<double>(subRows[0],subCols[0]));
    Buf.insert(Buf.end(),facDat[0]->begin(),facDat[0]->end());
    Buf.insert(Buf.end(),facDat[1]->begin(),facDat[1]->end());
    packBasic(Info,Buf);
}

void OperatorTensorProduct::addDat(const std::vector<std::complex<double> > &Buf, unsigned int Siz0, unsigned int Siz1){
    facDat.push_back(addComplex(hashString(subRows[0],subCols[0]),vector<complex<double> >(Buf.begin()+1,Buf.begin()+1+Siz0)));
    facDat.push_back(addComplex(hashString(subRows[1],subCols[1]),vector<complex<double> >(Buf.begin()+1+Siz0,Buf.begin()+1+Siz0+Siz1)));
}
