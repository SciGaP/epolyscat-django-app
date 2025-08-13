// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef EIGENINVERSEITER_H
#define EIGENINVERSEITER_H

#include "eigenSolverAbstract.h"
#include <complex>

/// brute force inverse iteration (as Arpack is letting me down)
class EigenInverseIter : public EigenSolverAbstract
{
    double _epsilon; // tolerance for convergence criterion
    int _updates;    // how many updates
    int _npower;     // which power to apply
    void _compute();
    bool converged(const std::vector<std::complex<double>> Eval);
public:
    EigenInverseIter(double Epsilon=1.e-12 /** maximal error in last few iterations */,
                     int Updates=4 /** maximal number of resolvent updates */,
                     int Power=5   /** maximal number of applies w/o update */);
    EigenInverseIter& withGuesseigenvalue(std::complex<double> Eguess);
};

#endif // EIGENINVERSEITER_H
