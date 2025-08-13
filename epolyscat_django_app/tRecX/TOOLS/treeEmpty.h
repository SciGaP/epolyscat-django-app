// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef TREEEMPTY_H
#define TREEEMPTY_H

#include "tree.h"
#include "treeDerived.h"
#include "coefficients.h"

///< Emtpy tree w/o any data (for debug purposes)
class TreeEmpty: public TreeDerived<TreeEmpty>
{
public:
    TreeEmpty(const Coefficients * T)
    {
        for(int k=0;k<T->childSize();k++)
            childAdd(new TreeEmpty(T->child(k)));
    }
};

#endif // TREEEMPTY_H
