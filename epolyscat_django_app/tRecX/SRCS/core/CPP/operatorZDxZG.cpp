// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorZDxZG.h"

#include "tools.h"
#include "useMatrix.h" // this includes the Eigen, which the compiler otherwise does not find

using namespace std;
#include "eigenNames.h"


OperatorZDxZG::OperatorZDxZG(std::vector<const UseMatrix*> Dat, string Kind)
    :OperatorTensorProduct(Dat,"ZDxZG"){

    string hash;

    vector<complex<double> > dDat(min(Dat[0]->rows(),Dat[0]->cols()));
    for(unsigned int d=0;d<dDat.size();d++)dDat[d]=(*Dat[0])(d,d).complex();

    hash=Kind+tools::str(Dat[0]->rows())+tools::str(Dat[0]->cols());
    facDat.push_back(addComplex(hash,dDat));

    vector<complex<double> > mDat(Dat[1]->rows()*Dat[1]->cols());
    unsigned int k=0;
    for(unsigned int n=0;n<Dat[1]->cols();n++)
        for(unsigned int m=0;m<Dat[1]->rows();m++,k++)
            mDat[k]=Dat[1]->operator()(m,n).complex();

    hash=Kind+tools::str(Dat[1]->rows())+tools::str(Dat[1]->cols());
    facDat.push_back(addComplex(hash,mDat));
}

//void OperatorZDxZG::axpy(std::complex<double> Alfa, const std::vector<std::complex<double> > &X,
//                         std::complex<double> Beta, std::vector<std::complex<double> > &Y) const {

//    // after some consulting of Eigen advice, Y +=Alfa * Matrix * X seems to be the right thing to do
//    if(Alfa==0.){
//        scale(Beta,Y);
//        return;
//    }

//    if(&X==&Y)ABORT("for now, no aliasing allowed");

//    vector<complex<double> >* cX=const_cast<vector<complex<double> >*>(&X);
//    if(Beta!=1.)scale(Beta,Y);

//    // most general case first
//    // NOTE: eigen native storage is colmn major; as our data is row major, left eigen index corresponds to right tRecX index
//    // NOTE: for efficiency, it may be worth while to choose a good sequence of applying 0 and 1
//    // NOTE: for efficiency, it may be useful to transpose one of the factors (check)
//    Map<MatrixXcd>(Y.data(),subRows[1],subRows[0]).noalias()
//            +=Alfa
//            * Map<MatrixXcd>(dat[1]->data(),subRows[1],subCols[1])
//            * Map<const MatrixXcd>(cX->data(),subCols[1],subCols[0])
//            * (Map<VectorXcd>(dat[0]->data(),dat[0]->size())).asDiagonal();
//}
void OperatorZDxZG::axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
                         const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const{

    // after some consulting of Eigen advice, Y +=Alfa * Matrix * X seems to be the right thing to do
    if(Alfa==0.){
        scale(Beta,Y,SizY);
        return;
    }

    if(X==Y)ABORT("for now, no aliasing allowed");

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
            * (Map<VectorXcd>(facDat[0]->data(),facDat[0]->size())).asDiagonal();
}

