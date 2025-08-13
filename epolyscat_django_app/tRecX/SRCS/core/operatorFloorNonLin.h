// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORFLOORNONLIN_H
#define OPERATORFLOORNONLIN_H

#include "operatorFloor.h"

class OperatorFloorNonLin: public OperatorFloor
{
public:
    OperatorFloorNonLin(std::string Kind):OperatorFloor(Kind){}
    virtual void updateNonLin(double Time, const Coefficients* C)=0;
};

#endif // OPERATORFLOORNONLIN_H
