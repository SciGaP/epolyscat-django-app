// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORZGXZG_H
#define OPERATORZGXZG_H


#include <string>
#include <vector>
#include <complex>

class UseMatrix;
class Index;

#include "qtEigenDense.h"
#include "operatorTensorProduct.h"

/** \ingroup OperatorFloors */
/// tensor product of complex general matrices
class OperatorZGxZG: public OperatorTensorProduct
{
public:
    OperatorZGxZG(std::vector<const UseMatrix *> Mat, std::string Kind);
    void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
                      const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const;

    OperatorZGxZG(const std::vector<int> &Info, const std::vector<std::complex<double> >&Buf)
        :OperatorTensorProduct(Info,Buf,"ZGxZG"){addDat(Buf,subRows[0]*subCols[0],subRows[1]*subCols[1]);}

    Eigen::MatrixXcd matrixFactor(int D) const;
};


#endif // OPERATORZGXZG_H
