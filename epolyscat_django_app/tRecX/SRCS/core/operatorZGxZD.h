// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORZGXZD_H
#define OPERATORZGXZD_H

#include <string>
#include <vector>
#include <complex>
#include "qtEigenDense.h"

class UseMatrix;
class Index;

#include "operatorTensorProduct.h"

/** \ingroup OperatorFloors */
/// tensor product of (complex diagon) x (complex general) matrices
class OperatorZGxZD: public OperatorTensorProduct
{
public:
    OperatorZGxZD(std::vector<const UseMatrix *> Mat, std::string Kind);
    void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
              const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const;
    OperatorZGxZD(const std::vector<int> &Info, const std::vector<std::complex<double> >&Buf)
        :OperatorTensorProduct(Info,Buf,"ZGxZD"){addDat(Buf,subRows[0]*subCols[0],subRows[1]);}
  Eigen::MatrixXcd matrixFactor(int D) const;
};

#endif // OPERATORZGXZD_H
