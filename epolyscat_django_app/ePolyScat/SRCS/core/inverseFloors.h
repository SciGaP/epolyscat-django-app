// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INVERSEFLOORS_H
#define INVERSEFLOORS_H

#include "inverse.h"
#include "operatorTree.h"

/// local inverse of OperatorTree: every diagonal block of Ovr is inverted
class InverseFloors : public Inverse
{
public:
    InverseFloors(const OperatorTree* Ovr);
    virtual void applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
    {Y.makeContinuous();};
    virtual void apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
    {OperatorTree::apply(A,Vec,B,Y);}

    /// for parallel application: use CoefficientsLocal
    virtual void applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const
    {Y.makeContinuous();};
    virtual void apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const
    {OperatorTree::apply(A,Vec,B,Y);}
};

#endif // INVERSEFLOORS_H
