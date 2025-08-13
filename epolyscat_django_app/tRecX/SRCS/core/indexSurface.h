// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXSURFACE_H
#define INDEXSURFACE_H

#include "indexDerived.h"

/// \ingroup
/// \brief special index for tSurff surfaces
class IndexSurface : public IndexDerived
{
public:
    IndexSurface(const Index* I, std::vector<double> Radius, unsigned int NSurf);
};

#endif // INDEXSURFACE_H
