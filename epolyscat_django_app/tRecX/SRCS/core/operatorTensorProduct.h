// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORTENSORPRODUCT_H
#define OPERATORTENSORPRODUCT_H

#include <string>
#include <vector>
#include <complex>

class UseMatrix;

#include "operatorFloor.h"

/// abstract base for all tensor products
class OperatorTensorProduct: public OperatorFloor{
protected:
    std::vector<std::shared_ptr<std::vector<std::complex<double>>>>facDat; ///< pointers for tensor factors
    std::vector<unsigned int> subRows,subCols; // dimension information of tensor factors
public:
    OperatorTensorProduct(std::vector<const UseMatrix *> Mat, std::string Kind);
    OperatorTensorProduct(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf, std::string Kind);
    void pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const;
    void addDat(const std::vector<std::complex<double> > &Buf, unsigned int Siz0,unsigned int Siz1);
};


#endif // OPERATORTENSORPRODUCT_H
