// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef FINDIFF_H
#define FINDIFF_H

#include <complex>
#include "useMatrix.h"

class Algebra;

/// \ingroup Basissets
/// @brief finite difference schemes (not currently maintained)
///
/// set up finite difference schemes for a choice of operators and interpolation assumptions
class FinDiff
{

    std::vector<double> grid,weig;

    // scaling functions, parameters, and internal state
    Algebra *argScal,*derScal,*scalD1,*scalD2;
    double r0;
    std::complex<double> eta;
    std::vector<double> singularX;
    double xLast;
    std::complex<double> scalLast;

    // scaling function and its derivative
    std::complex<double> scalArg(double X);
    std::complex<double> scalDer(double X);

    // functions to base the finite difference scheme on (usually: polynomials)
    void interFunc(unsigned int Points, std::complex<double> Z, std::complex<double> A, std::complex<double> B,
                          std::vector<std::vector<std::complex<double> > >& V);

    void construct(const std::string Kind="smooth");
    void scheme(unsigned int Points, UseMatrix & Der, UseMatrix &Lap); ///< 1d FD schemes for various operators

public:
    static std::string FDgrid; ///< external control for FD grid type
    static std::string FDmatrix;  ///< how to construct FD matrix

    FinDiff(){}

    /// arbitrary input grid including weights, if nonequidistant
    FinDiff(UseMatrix Grid, UseMatrix Weight=UseMatrix(), double R0=0., std::complex<double> Eta=1.,const std::string Kind="sudden");

    /// standard equidistant grid, NOT containing end-points X0,X1
    FinDiff(double X0, double X1, unsigned int N, double R0=0., std::complex<double> Eta=1., const std::string Kind="sudden");

    /// finite differences matrices for various operators
    static void matrix(std::string Oper, unsigned int Points, const UseMatrix & Grid, const UseMatrix & Weig,
                       double R0, std::complex<double> Eta, double Par, UseMatrix & Matrix, std::string Kind);
    static void test();
};

#endif // FINDIFF_H
