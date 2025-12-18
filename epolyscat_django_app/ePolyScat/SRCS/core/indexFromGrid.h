// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXFROMGRID_H
#define INDEXFROMGRID_H

#include "index.h"

/// construct Index by converting grids to basis where desired and possible
///
/// default bases for given coordinate are used (where no basis kind is specified explicitly)
class IndexFromGrid: public Index
{
    static const BasisAbstract * basisFromGrid(const Index *GridIdx, std::string ContractFactor, std::vector<const Index *> Path);
    IndexFromGrid(const Index* IdxGrid, std::vector<std::string> Deflate, std::vector<const Index *> Path);
public:
    IndexFromGrid(const Index* IdxGrid, std::vector<std::string> Deflate={} /** reduce basis size from grid size by Deflate  at given level */)
        :IndexFromGrid(IdxGrid,Deflate,{}){}
};

#endif // INDEXFROMGRID_H
