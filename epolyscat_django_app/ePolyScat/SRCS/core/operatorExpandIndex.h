// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATOREXPANDINDEX_H
#define OPERATOREXPANDINDEX_H

#include "operatorTree.h"

/// OperatorTree from Factor and Identity on tensor product space
/// allowing constraints among factor indices
///
/// Algorithm:
/// - descends Factor
/// - if MatchIdx matches Idx: branches as in Factor
/// - else                   : branches for Id on Idx level
/// a new index is on the other side of Idx is created
/// for correct function, all indices on Match of Factor must appear in Idx
class OperatorExpandIndex: public OperatorTree
{
    static void coefficients(const Coefficients* Fac, Coefficients *C);
public:
    /// expands Factor for idx() (MatchI=true) or else jdx() to match Mdx
    /// create new index on the other side
    OperatorExpandIndex(const OperatorTree* Factor, const Index* Mdx, bool MatchI /** =true for lhs(=I),  else  rhs(=J) */);

    /// tensor expand coefficients by 1's into Coefficients(Idx)
    static Coefficients coefficients(const Coefficients& Fac, const Index* IdxX);
    static bool debug;
};

#endif // OPERATOREXPANDINDEX_H
