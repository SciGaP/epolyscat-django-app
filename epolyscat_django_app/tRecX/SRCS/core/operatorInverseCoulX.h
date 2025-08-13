// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORINVERSECOULX_H
#define OPERATORINVERSECOULX_H

#include "operatorTree.h"

class Index;

class OperatorInverseCoulX:public OperatorTree
{
public:
    OperatorInverseCoulX(const Index *Idx, unsigned int SubD, unsigned int SuperD, bool BandOvr);
};

#endif // OPERATORINVERSECOULX_H
