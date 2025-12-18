// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORZDXZD_H
#define OPERATORZDXZD_H


#include <string>
#include <vector>
#include <complex>

class UseMatrix;
class Index;
class OpeatorZD;

#include "qtEigenDense.h"
#include "operatorTensorProduct.h"
#include "operatorZD.h"

/** \ingroup OperatorFloors */
/// tensor product of diagonal matrices (returns factors separately)
class OperatorZDxZD: public OperatorZD
{
    std::vector<OperatorZD> facs;
public:
    OperatorZDxZD(std::vector<const UseMatrix *> Mat, std::string Kind);
    Eigen::MatrixXcd matrixFactor(int D) const;
};


#endif // OPERATORZDXZD_H
