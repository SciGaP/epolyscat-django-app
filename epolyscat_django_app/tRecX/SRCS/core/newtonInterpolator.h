// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef NEWTONINTERPOLATOR_H
#define NEWTONINTERPOLATOR_H
#include "toolsHeader.h"

class NewtonInterpolator {
  public:
    NewtonInterpolator(std::vector<std::complex<double> >& coordinates, std::vector<std::complex<double> >& fctValues);
    std::vector<std::complex<double> >& getVals(std::vector<std::complex<double> >& coordinates, std::vector<std::complex<double> >& fctValues) const;
  protected:
    std::vector<std::complex<double> > coeffs;
    std::vector<std::complex<double> > supports;
};

class Wavefunction;
class Discretization;
class Index;

class NewtonInterpolatorWF {
    // dividierte differenzen schema
    // knownWfs->time represent the support points
    // knownWfs->coefs represent the known function values at the support points
    // the interpolation at a given time returns a wavefunction
public:
    NewtonInterpolatorWF(const Index* Idx, int numberOfSupportPoints); // allocates some memory
    NewtonInterpolatorWF(const Discretization* disc, int numberOfSupportPoints); // allocates some memory
    ~NewtonInterpolatorWF();

    void getInterpolatedWF(double time, Wavefunction* result); // computePolynomialCoefficients needs to be called first ! computes the interpolated wavefunction
    void timeDerivative(double time, Wavefunction* result); // computePolynomialCoefficients needs to be called first ! computes the interpolated wavefunction

    void computePolynomialCoefficients(const std::vector<Wavefunction*> & knownWfs); // for given support points compute the expansion coefficients of the unique polynomial in newton basis
    bool inInterval(double Time);

    std::vector<Wavefunction*> polynomialCoefficients; // the expansion coefficients of the polynomial
    std::vector<Wavefunction*> divDiffs; // helper variable for computePolynomialCoefficients
    std::vector<double> times; // the supportpoints
};

#endif // NEWTONINTERPOLATOR_H
