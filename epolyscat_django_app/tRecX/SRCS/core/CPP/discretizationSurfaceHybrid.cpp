// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "discretizationSurfaceHybrid.h"

#include "discretizationHybrid.h"
#include "operatorTree.h"
#include "basicDisc.h"
#include "printOutput.h"

using namespace std;

DiscretizationSurfaceHybrid::DiscretizationSurfaceHybrid(const DiscretizationHybrid *Parent, const std::vector<double> &Rad, unsigned int NSurf)
    :DiscretizationSurface(Parent->comp[1],Rad,NSurf)
{
    if(idx()->axisName()!="Neut&Chan" and idx()->axisName()!="Subspace&Complement")
        PrintOutput::warning("not tested for hybrid discretization "+idx()->axisName());

    // only adjust name and map
    _mapFromParent.reset(new Map(mapFromParent(),Parent->idx()));
    name="(Hybrid)"+name;
}
