// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "inverseFloors.h"
#include "operatorFloor.h"

InverseFloors::InverseFloors(const OperatorTree* Ovr)
    :Inverse("Inv("+Ovr->name+")",Ovr->iIndex,Ovr->jIndex)
{
    OperatorTree::nodeCopy(Ovr,false);
    if(OperatorTree::iIndex==OperatorTree::jIndex and floor())floorInvert();
    for(int k=0;k<Ovr->childSize();k++)Inverse::childAdd(new InverseFloors(Ovr->child(k)));
}
