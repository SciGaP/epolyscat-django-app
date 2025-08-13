// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INVERSEDVRMIXFE_H
#define INVERSEDVRMIXFE_H

#include "inverse.h"
#include "tree.h"

class Index;
class OperatorTree;

/// inverse of a DVR basis that is interspersed with isolated FE (full) overlap matrices
///
/// (this appears to be working)
class InverseDVRmixFE: public Inverse
{
    OperatorTree* inv0;
public:
    InverseDVRmixFE(Index *Idx);
    void apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    void applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

    void applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const;
    void apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const;
    void parallelSetup() const;// {ABORT("no parallel setup implmented");}
};

#endif // INVERSEDVRMIXFE_H
