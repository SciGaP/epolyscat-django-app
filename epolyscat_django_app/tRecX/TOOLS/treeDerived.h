// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef TREEDERIVED_H
#define TREEDERIVED_H

#include "tree.h"

template <typename T>
class TreeDerived: public Tree<T>{
public:
    TreeDerived(const T* Parent=0):Tree<T>(Parent){} ///< elementary constructor
    TreeDerived(const TreeDerived<T> &Other, bool View=false); ///< (deep) copy constructor
    TreeDerived<T> & operator=(const TreeDerived<T> & Other); ///< (deep) assignement operator
};
#endif // TREEDERIVED_H
