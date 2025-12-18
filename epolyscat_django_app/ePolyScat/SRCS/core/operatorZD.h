// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORZD_H
#define OPERATORZD_H


#include <string>
#include <vector>
#include <complex>

class UseMatrix;
class Index;

#include "operatorFloor.h"

/** \ingroup OperatorFloors */
/// complex diagonal matrix
class OperatorZD: public OperatorFloor
{
protected:
    void axpy(const std::complex<double> & Alfa, const std::complex<double>*X, unsigned int SizX,
              const std::complex<double> & Beta,       std::complex<double>*Y, unsigned int SizY) const;

    void axpyTranspose(std::complex<double> Alfa, const std::vector<std::complex<double> > & X,
                       std::complex<double> Beta, std::vector<std::complex<double> > & Y) const
    {axpy(Alfa,X.data(),X.size(),Beta,Y.data(),Y.size());}
    void construct(const UseMatrix *Mat, std::string Kind);
    void construct(const Eigen::VectorXcd & Mat, std::string Kind);
public:
    OperatorZD():OperatorFloor(0,0,"ZD"){}
    OperatorZD(const UseMatrix *Mat, std::string Kind);
    OperatorZD(const Eigen::MatrixXcd & Mat, std::string Kind);
    OperatorZD(const std::complex<double>* Diag, int Size);

    void pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const;

    OperatorZD(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf):OperatorFloor("ZD"){
        unpackBasic(Info,Buf);
        dat=addComplex(hashString(_rows,_cols),std::vector<std::complex<double> >(Buf.begin(),Buf.begin()+Info[3]));
    }
    
    virtual long applyCount() const;
};


#endif // OPERATORZD_H
