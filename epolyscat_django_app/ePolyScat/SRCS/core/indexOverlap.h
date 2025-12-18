// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXOVERLAP_H
#define INDEXOVERLAP_H

#include <vector>

class Index;
class OperatorAbstract;
class OperatorTree;
class Inverse;

/// manage the overlap and inverse for Index
class IndexOverlap
{
    IndexOverlap();
public:
    static void localOverlapAndInverse(const Index* Idx, OperatorTree* Ovr, OperatorTree* Inv);
    static const OperatorTree* getOverlap(const Index *Idx, std::vector<unsigned int> &idx, const Index *&iRoot, int Kind);
    static void set(const Index *Idx, const OperatorAbstract *Ov=0);
    static const OperatorAbstract* get(const Index* Idx);
    static OperatorTree * getTree(const Index* Idx, int Kind);
    static void setInverseOverlap(const Index* Idx,const Inverse *Inv);
    static const Inverse * inverseOverlap(const Index* Idx);
};
#endif // INDEXOVERLAP_H
