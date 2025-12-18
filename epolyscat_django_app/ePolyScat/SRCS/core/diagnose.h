// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DIAGNOSE_H
#define DIAGNOSE_H

#include <climits>

class OperatorAbstract;

class Diagnose
{
public:
    Diagnose();
    void show(const OperatorAbstract* Op, unsigned int BlockLevel=INT_MAX) const;
};

#endif // DIAGNOSE_H
