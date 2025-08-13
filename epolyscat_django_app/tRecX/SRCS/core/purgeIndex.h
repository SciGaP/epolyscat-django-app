// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PURGE_INDEX_H
#define PURGE_INDEX_H

#include <vector>
#include <set>
#include <utility>

class Index;
class OperatorTree;

/**
 * Remove unused children from index structure. This need arises mainly when working with constrained
 * indices. An OperatorTree from a constrained index to a newly generated one might not map into all
 * subspaces. Given an index and one or many OperatorTrees, this class removes all children from the index,
 * that are not mapped into/from by any OperatorTree.
 */
class PurgeIndex{
    Index* idx;

    // the bool value indicates: true <=> iIndex, false <=> jIndex
    std::vector<std::pair<const OperatorTree*, bool>> opTrees;

    std::set<const Index*> used;

    void registerUsed(const OperatorTree* Op, bool IIndex);
    bool purge(Index* Idx);

public:
    PurgeIndex(Index* Idx);
    PurgeIndex& usingOperatorTree(const OperatorTree* Op);
    void run();

    /// prepare before run is optional
    void prepare();

    /// isUsed only after prepare
    bool isUsed(const Index* Idx) const;
};



#endif
