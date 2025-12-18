// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef VOLKOV_GRID_H
#define VOLKOV_GRID_H

#include "operatorAbstract.h"
#include "discretizationDerived.h" //temporary

class Discretization;
class DiscretizationSurface;
class DiscretizationDerived;
class SurfaceFlux;
class CoefficientsFunction;

class VolkovGrid: public OperatorAbstract{
    DiscretizationDerived* gridSpecDisc;
    CoefficientsFunction* volkov;
    Coefficients* cGrid;
    double _time;

public:
    VolkovGrid(const Discretization * D, SurfaceFlux *Flux, const DiscretizationDerived* SmoothSpecDisc):VolkovGrid(Flux,SmoothSpecDisc->idx()){}
    VolkovGrid(SurfaceFlux *Flux, const Index *SpecIdx);

    void update(double Time, const Coefficients* CurrentVec);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const override;
};


#endif
