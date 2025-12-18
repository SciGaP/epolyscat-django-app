// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORVECTORS_H
#define OPERATORVECTORS_H

#include <vector>
#include <string>

#include "operatorAbstract.h"
#include "coefficients.h"

class Coefficients;
class Index;
/** \ingroup Structures */
/// set of vectors considered as columns of an operator
class OperatorVectors : public OperatorAbstract
{
    std::vector<Coefficients> vecs;
public:
    OperatorVectors(std::string Name, const Index* IIndex, const Index* JIndex);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    void insertColumn(unsigned int J, const Coefficients & C);
};

#endif // OPERATORVECTORS_H
