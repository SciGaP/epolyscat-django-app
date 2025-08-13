// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORWITHINVERSE_H
#define OPERATORWITHINVERSE_H

#include "operatorAbstract.h"
#include "index.h"
#include "inverse.h"
#include "coefficients.h"

/// maps to Parent space (assuming original operator maps to Dual)
class OperatorWithInverse : public OperatorAbstract
{
    std::shared_ptr<const OperatorAbstract> _op;
    mutable Coefficients tmp;
public:
    OperatorWithInverse(std::shared_ptr<const OperatorAbstract> Op)
        :OperatorAbstract(Op->name+"_withInverse",Op->iIndex,Op->jIndex)
    {
        _op=Op;
        tmp.reset(_op->iIndex);
    }

    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
        _op->apply(A,Vec,0.,tmp);
        if(not iIndex->inverseOverlap())DEVABORT("iIndex does not have inverse");
        iIndex->inverseOverlap()->apply(1.,tmp,B,Y);
    }
};

#endif // OPERATORWITHINVERSE_H
