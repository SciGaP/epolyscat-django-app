// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "convertFromGrid.h"

#include "index.h"
#include "operatorMap.h"

using namespace std;

ConvertFromGrid::ConvertFromGrid(const Index *GridIdx, std::vector<int> Deflate)
    :gridIdx(GridIdx)
{
    if(Deflate.size()==0)Deflate.assign(GridIdx->height(),1);
    basIdx=gridIdx->toIndexBasis(Deflate);
    toBas=new OperatorMap(basIdx,gridIdx);
}
