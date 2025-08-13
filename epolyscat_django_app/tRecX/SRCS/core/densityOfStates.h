// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DENSITYOFSTATES_H
#define DENSITYOFSTATES_H

#include <vector>
#include <complex>
#include <string>

class OperatorAbstract;
class Coefficients;
class ReadInput;

class DensityOfStates
{
    double eMin,eMax;
    int nPoints;
    double box;
    std::string outDir;
public:
    DensityOfStates(ReadInput & Inp);
    void output(const OperatorAbstract* Op);
    std::vector<double> compute(const OperatorAbstract* Op);
};

#endif // DENSITYOFSTATES_H
