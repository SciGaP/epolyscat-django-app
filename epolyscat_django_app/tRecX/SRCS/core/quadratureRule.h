// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef QUADRATURERULE_H
#define QUADRATURERULE_H

#include <vector>
#include <complex>

class Algebra;
class BasisIntegrable;
///@brief quadratur rule for BasisIntegrable
namespace QuadratureRule
{

void pointsAndWeights(const BasisIntegrable * Bas, std::vector<double> & Points, std::vector<double> & Weights);

std::vector<std::complex<double> > integralsBasisAlgebra(const BasisIntegrable * Bas, const Algebra* Alg);
std::vector<std::complex<double> > integralsBasisBasis(const BasisIntegrable * Bas);

void test();
}

#endif // QUADRATURERULE_H
