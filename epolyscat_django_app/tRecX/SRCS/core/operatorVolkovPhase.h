// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORVOLKOVPHASE_H
#define OPERATORVOLKOVPHASE_H

#include "operatorAbstract.h"

class DiscretizationDerived;
class CoefficientsFunction;

class OperatorVolkovPhase:public OperatorAbstract
{
    CoefficientsFunction* volkov;
    void update(double Time, const Coefficients* CurrentVec);

public:
    OperatorVolkovPhase(const DiscretizationDerived *GridSpecDisc, const std::string KLevel);
    void apply(const std::complex<double> A, const Coefficients& X, const std::complex<double> B, Coefficients& Y) const;
};

#endif // OPERATORVOLKOVPHASE_H
