// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef SUBITER_H
#define SUBITER_H
#include "tools.h"

class Coefficients;

/// @brief subspace iteration class for Eigenvalues
///
/// get eigenvectors and values of and Operator in the vicinity of a complex value E0
/// use an invertible Op0 to generate new vectors: newVectors = (Op0-E0)^-1 Operator oldVectors
class SubIter
{
public:
    SubIter();

    void eigen(unsigned int Nvec, std::vector<Coefficients*> Rvec, std::complex<double> Eval, double Eps=1.e-12);
};

#endif // SUBITER_H
