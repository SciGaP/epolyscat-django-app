// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMAT2D_H
#define BASISMAT2D_H

#include <string>
#include "qtEigenDense.h"
#include "index.h"

class BasisAbstract;
class UseMatrix;
class BasisMat2D
{
    Eigen::MatrixXcd _mat,_mat0,_mat1;
public:
    BasisMat2D(){}
    BasisMat2D(std::string Op, const Index * IIndex, const Index* JIndex);
    bool isEmpty() const {return _mat.size()==0 and (_mat0.size()==0 or _mat1.size()==0);}
//    const UseMatrix useMat() const;
    const std::vector<const Eigen::MatrixXcd*> mats() const;
};
#endif // BASISMAT1D_H
