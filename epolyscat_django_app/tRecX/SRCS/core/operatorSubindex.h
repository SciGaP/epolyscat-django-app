#ifndef OPERATORSUBINDEX_H
#define OPERATORSUBINDEX_H

#include "operatorTree.h"

class OperatorSubindex : public OperatorTree
{
public:
    /// Operator with only those blocks of Op that match SubIdx on the MatchI side,
    /// remove blocks not in SubIdx, replace index with SubIdx
    OperatorSubindex(const OperatorTree *Op, const Index* SubIdx, bool MatchI /** true: match idx, else: jdx */);
};

#endif // OPERATORSUBINDEX_H
