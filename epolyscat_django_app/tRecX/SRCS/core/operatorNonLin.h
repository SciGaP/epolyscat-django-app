// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORNONLIN_H
#define OPERATORNONLIN_H


#include "eigenSolverAbstract.h"
#include "eigenSolver.h"

class OperatorNonLin : public OperatorTree
{
protected:
    OperatorFloor * oFloor;
    bool _view;
public:
    OperatorNonLin(const std::string Name, const std::string & Definition, const Index *IIndex, const Index *JIndex);
    void update(double Time, Coefficients* C);
};
#endif
