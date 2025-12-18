// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "selectionRule.h"

#include "discretization.h"
#include "index.h"
#include "axis.h"

using namespace std;

// a rule name for coordinates a,b,c has format name:a.b.c
//map<string,string> SelectionRule::operList["XY:Phi.Eta"]="<sin(Q)*cos(Q)><1-Q*Q>";

void SelectionRule::add(const std::vector<Axis> &Ax, std::string Rule)
{
}


