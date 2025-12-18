// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXPROD_H
#define INDEXPROD_H

#include "index.h"

/// \ingroup Index
/// \brief tensor product of two indices
class IndexProd : public Index
{
public:
    IndexProd(const Index *A, const Index * B);
};

#endif // INDEXPROD_H
