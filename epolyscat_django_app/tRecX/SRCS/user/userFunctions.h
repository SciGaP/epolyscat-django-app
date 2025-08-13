// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef USERFUNCTIONS_H
#define USERFUNCTIONS_H

#include "basisMat.h"

/** @brief example user-defined function
 *
 * add new functions (if needed ) by this example
 * - define the class here
 * - add code to userFunctions.cpp file \n
 *   give name that will be used for input \n
 *   e.g. basisMatFunc::name="myFunc" can be addressed in input as
 *
 *   <myFunc(Q)> for single argument functions \n
 *   with Q a place-holder for current coordinate
 *
 *   <{}>...<myFunc{Phi,Y}> \n
 *   for a two-argument function taking the coordinates from axes Phi and Y;\n
 *   the place-holder <{}> must be in the position of one of the axes, say Phi\n
 *   the last factor must at the position of other axis, say Y \n
 *   the order Phi, Y corresponds to arguments X[0], X[1]
 *
 * - include the new class in UserFunctions::list()
 *
 *   # CAUTION:
 *   for now, names are not checked for uniqueness
 *   and multiple names can cause undefined behavior
 */


class basisMatUserRedefined:public basisMatFunc {
    int kind;
public:
    basisMatUserRedefined();
    std::complex<double> operator()(std::vector<std::complex<double> >X) const;
};

/// @brief example PyrA function
class basisMatPyrA:public basisMatFunc {
    int kind;
public:
    basisMatPyrA(int Kind);
    std::complex<double> operator()(std::vector<std::complex<double> >X) const;
};

/// @brief 2d potential test: harmonic x Morse
class basisMatHarMorse:public basisMatFunc {
    int kind;
    AlgebraMorse *morse0,*morse1;
public:
    ~basisMatHarMorse();
    basisMatHarMorse(int Kind);
    std::complex<double> operator()(std::vector<std::complex<double> >X) const;
};

/// @brief electron-electron repulsion 2X1d helium
class basisMat2Dhelium:public basisMatFunc {
public:
    basisMat2Dhelium();
    std::complex<double> operator()(std::vector<std::complex<double> >X) const;
};
/// @brief for debugging
class basisMatXtimesY:public basisMatFunc {
public:
    basisMatXtimesY();
    std::complex<double> operator()(std::vector<std::complex<double> >X) const;
};

#endif // USERFUNCTIONS_H
