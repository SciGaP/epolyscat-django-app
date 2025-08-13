// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef CONVERTFROMGRID_H
#define CONVERTFROMGRID_H

#include <vector>

class Index;
class OperatorMap;

class ConvertFromGrid
{
    const Index *gridIdx, *basIdx;
    const OperatorMap* toBas;
public:
    ConvertFromGrid(const Index * GridIdx, std::vector<int> Deflate={} /** reduce basis size at level k by Deflate[k] compared to grid size (if specified) */);
};

#endif // CONVERTFROMGRID_H
