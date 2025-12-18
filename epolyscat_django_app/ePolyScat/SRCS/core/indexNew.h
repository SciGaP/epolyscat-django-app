// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXNEW_H
#define INDEXNEW_H

#include "index.h"

class AxisTree;

///@brief temporary for developing new index constructor
class IndexNew : public Index
{
    /// add child, IF the subtree is non-empty
    void addNonEmpty(AxisTree *Ax, const IndexConstraint* Constraint, const std::vector<unsigned int> Pos,
                     const std::vector<const Index *> Path, std::vector<int> & Removed);
    void setBasisRemoved(const std::vector<int> & Removed); ///< handle basis truncation
    IndexNew(const AxisTree *Ax, const IndexConstraint* Constraint,
          std::vector<unsigned int> InPos,
          std::vector<const Index*> InPath, bool Entry);
public:
    static void read(ReadInput & Inp);

    virtual ~IndexNew();
    IndexNew():Index(){}
    IndexNew(const Index & Other):Index(Other){}

    static bool doNotUse;
    /// new basic construct from axes
    IndexNew(const AxisTree *Ax, const IndexConstraint* Constraint=0)
        :IndexNew(Ax,Constraint,{},{},true){}

    ///@brief interpretation of basis strings that depends on position in tree
    static void resolveDependence(BasisSetDef & Def,std::vector<unsigned int> Branch, const std::vector<const Index*> Path);

    ///@brief set up overlap and inverse for Index's
    void buildOverlap();
    const AxisTree* axes() const;
};



#endif // INDEXNEW_H
