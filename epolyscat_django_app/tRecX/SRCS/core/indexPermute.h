// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXPERMUTE_H
#define INDEXPERMUTE_H

#include "index.h"

class IndexPermute : public Index
{
    std::vector<size_t> _permC;
public:
    IndexPermute(const Index* Idx, std::string Hierarchy);
    const std::vector<size_t> & permC() const { return _permC;};
};

#endif // INDEXPERMUTE_H
