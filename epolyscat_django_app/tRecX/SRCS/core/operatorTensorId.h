// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORTENSORID_H
#define OPERATORTENSORID_H

#include "operatorAbstract.h"
#include <vector>

class IndexProd;
class Coefficients;
class CoefficientsPermute;


/// Op = Id (x) OpFactor, permuted to match rhs index JFull
class OperatorTensorId : public OperatorAbstract
{
    const OperatorAbstract* factor;
    std::vector<Coefficients*> iVfac,jVfac; // points to factor level
    CoefficientsPermute *iLoc,*jLoc; // local storage with product structure; provides permuted view matching original structure
    IndexProd *iProd,*jProd; // product indices
    void indexReplace(Coefficients* C, const Index *I);
public:
    OperatorTensorId(const OperatorAbstract * Factor, const Index *JFull, const Index *IFull=0);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
};

#endif // OPERATORTENSORID_H
