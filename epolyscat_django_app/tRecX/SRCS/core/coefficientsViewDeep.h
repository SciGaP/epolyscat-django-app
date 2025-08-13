// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COEFFICIENTSVIEWDEEP_H
#define COEFFICIENTSVIEWDEEP_H

#include <string>
#include <complex>

#include "coefficients.h"

class CoefficientsViewDeep {
    Coefficients _view;
    static void extend(Coefficients* View, std::complex<double>* Data);
    static void view(Coefficients* View, Coefficients* C);
    static void disown(Coefficients*);
public:
    CoefficientsViewDeep(const Index* Idx, int Depth=INT_MAX, bool paraSub=false);
    Coefficients* view(Coefficients* C);
    std::string str() const {return _view.str();}
};

#endif // COEFFICIENTSVIEWDEEP_H
