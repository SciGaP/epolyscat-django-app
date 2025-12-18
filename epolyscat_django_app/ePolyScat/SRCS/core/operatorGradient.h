// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORGRADIENT_H
#define OPERATORGRADIENT_H

#include "operatorAbstract.h"

class Discretization;
class DiscretizationGrid;
class ReadInput;

/** \ingroup Structures */
/// values and partial derivatives on a grid
class OperatorGradient : public OperatorAbstract
{
    DiscretizationGrid * grid;
    Coefficients * model,*permModel;

public:
    /// return pointer to operator gradient from Parent with grid read from Inp
    static OperatorGradient* read(const Discretization * Parent /** original discretization */,
                                  ReadInput & Inp /** read grid specification from here */);
    OperatorGradient(const Discretization *Parent, std::vector<std::string> &Ax,
                     std::vector<unsigned int> &Points, std::vector<std::vector<double> > &Bounds);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
};

#endif // OPERATORGRADIENT_H
