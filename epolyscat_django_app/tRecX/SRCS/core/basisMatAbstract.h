// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMATABSTRACT_H
#define BASISMATABSTRACT_H

#include "qtEigenDense.h"
#include "index.h"

class BasisMatAbstract
{
protected:
    Eigen::MatrixXcd _mat;
    virtual void _construct(std::string Op, const Index * IIndex, const Index* JIndex)=0;
public:
    BasisMatAbstract(){}
    bool isEmpty() const {return _mat.size()==0;}
    const UseMatrix useMat() const;
    const std::vector<const Eigen::MatrixXcd*> mats() const;
public:
};

#endif // BASISMATABSTRACT_H
