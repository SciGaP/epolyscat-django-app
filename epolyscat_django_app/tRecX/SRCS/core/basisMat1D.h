// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMAT1D_H
#define BASISMAT1D_H

#include <string>
#define EIGEN_MATRIX_PLUGIN "EigenAddonMatrix.h"
#include "index.h"
#include "qtEigenDense.h"

#include "basisMatMatrix.h"

class BasisAbstract;
class UseMatrix;
/// matrix for 1d BasisIntegrable
class BasisMat1D : public BasisMatMatrix
{
    Eigen::MatrixXcd valDer(const UseMatrix & X, const BasisIntegrable* Bas, int Derivative) const;
public:
    BasisMat1D(){}
    BasisMat1D(std::string Op, int Idim, int Jdim);
    BasisMat1D(std::string Op, const Index * IIndex, const Index* JIndex);
    BasisMat1D(std::string Op, const BasisAbstract* IBas, const BasisAbstract * JBas);
};
#endif // BASISMAT1D_H
