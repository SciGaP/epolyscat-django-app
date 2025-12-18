// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXASSEMBLED_H
#define INDEXASSEMBLED_H

#include "index.h"

/// Index that relates to original Index
class IndexAssembled : public Index
{
    const Index* _from;
public:
    IndexAssembled(const Index* From):_from(From){nodeCopy(From,false);}
    const Index* from() const {return _from;};
};

#endif // INDEXASSEMBLED_H
