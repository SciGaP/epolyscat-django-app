// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorZGxZG.h"
#include "qtEigenDense.h"
#include "tools.h"
#include "useMatrix.h" // this includes the Eigen, which the compiler otherwise does not find

using namespace std;
#include "eigenNames.h"


OperatorZGxZG::OperatorZGxZG(std::vector<const UseMatrix*> Dat, string Kind)
    :OperatorTensorProduct(Dat,"ZGxZG"){

    string hash;

    for(unsigned int i=0;i<2;i++){
        vector<complex<double> > mDat(Dat[i]->rows()*Dat[i]->cols());

        unsigned int k=0;
        for(unsigned int n=0;n<Dat[i]->cols();n++)
            for(unsigned int m=0;m<Dat[i]->rows();m++,k++)
                mDat[k]=Dat[i]->operator()(m,n).complex();

        hash=Kind+tools::str(Dat[i]->rows())+tools::str(Dat[i]->cols());
        facDat.push_back(addComplex(hash,mDat));
    }
}

Eigen::MatrixXcd OperatorZGxZG::matrixFactor(int D) const{
    return Map<MatrixXcd>(facDat[D]->data(),subRows[D],subCols[D]);
}

void OperatorZGxZG::axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
                         const std::complex<double> & Beta,       std::complex<double>* Y, unsigned int SizY) const{
    // after some consulting of Eigen advice, Y +=Alfa * Matrix * X seems to be the right thing to do
    if(Alfa==0.){
        scale(Beta,Y,SizY);
        return;
    }

    if(&X==&Y)ABORT("for now, no aliasing allowed");

    complex<double>* cX=const_cast<complex<double>*>(X);
    if(Beta!=1.)scale(Beta,Y,SizY);

    // most general case first
    // NOTE: eigen native storage is colmn major; as our data is row major, left eigen index corresponds to right tRecX index
    // NOTE: for efficiency, it may be worth while to choose a good sequence of applying 0 and 1
    // NOTE: for efficiency, it may be useful to transpose one of the factors (check)
    Map<MatrixXcd>(Y,subRows[1],subRows[0]).noalias()
            +=Alfa
            * Map<MatrixXcd>(facDat[1]->data(),subRows[1],subCols[1])
            * Map<const MatrixXcd>(cX,subCols[1],subCols[0])
            * (Map<MatrixXcd>(facDat[0]->data(),subRows[0],subCols[0])).transpose();
}
