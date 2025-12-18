// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXQUOT_H
#define INDEXQUOT_H

#include "index.h"

/// \ingroup Index
/// \brief builds a quotient index - the levels NOT contained in Fac
class IndexQuot : public Index
{
public:
    IndexQuot(const Index * Full, const Index * Fac);
};

#endif // INDEXQUOT_H
