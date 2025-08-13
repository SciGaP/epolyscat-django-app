// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INVERSEDVR_H
#define INVERSEDVR_H

#include "coefficients.h"
#include "inverse.h"
#include "tree.h"

class CoefficientsLocal;
class Index;
class OperatorFloor3d;


class InverseDVR : public Inverse
{
    const Coefficients* _diagonal;
    const Coefficients* diagonalLocal;
    InverseDVR(const Coefficients * Diagonal);
public:
    ~InverseDVR(){delete _diagonal;}
    // note: here we need the Discretization explicitly as no Operator constructor is not optimal
    InverseDVR(const Index *Idx);

    // only make continuous
    void applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
    {Y.makeContinuous();}

    void apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

    /// create the diagonalLocal view on diagonal
    void parallelSetup() const;
    void apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const;
    void applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const
    {Y.makeContinuous();}

    std::complex<double> invWeig(int K) const {if(not isLeaf())ABORT("only on floor level for now");return const_cast<Coefficients*>(_diagonal)->data()[K];}

    const std::complex<double> * diagonal() const {return _diagonal->data();}
};


#endif // INVERSEDVR_H
