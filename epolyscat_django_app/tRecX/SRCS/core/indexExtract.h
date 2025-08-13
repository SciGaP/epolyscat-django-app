// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXEXTRACT_H
#define INDEXEXTRACT_H

#include "index.h"
typedef bool (*selectIndex)(const Index * Idx); /// slection criterion

/// \ingroup Index
/// \brief Extract new from Index by various criteria
class IndexExtract : public Index
{
    void _constructOld(Index* Result, const Index *Idx, selectIndex Select);
    Index* _construct(const Index *Idx, selectIndex Select);
    IndexExtract(const Index* Idx, const std::vector<std::string> Ax, bool Complement, const IndexConstraint* Constraint, bool Entry);
public:
    ///\brief all nodes with axisName in Ax (or complement), avoid double-count when extracting from product space, re-impose constraints, if any
    IndexExtract(const Index* Idx, const std::vector<std::string> Ax, bool Complement, const IndexConstraint* Constraint=0)
        :IndexExtract(Idx,Ax,Complement,Constraint,true){}
    ///\brief return new IndexExtract(...) or 0 if emtpy
    static Index* get(const Index* Idx, const std::vector<std::string> Ax, bool Complement, const IndexConstraint* Constraint=0);
};

#endif // INDEXEXTRACT_H
