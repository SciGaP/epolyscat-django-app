// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef AXISTREE_H
#define AXISTREE_H

#include "tree.h"
#include "axis.h"

class ReadInput;
class AxisTree : public Axis, public Tree<AxisTree>
{
    std::string _subset;
    void insertHybrid();
public:
    virtual ~AxisTree(){}
    AxisTree():_subset(""){}
    AxisTree(const Axis Ax):Axis(Ax),_subset(""){}
    AxisTree(const std::vector<Axis> & Axes);
    AxisTree(ReadInput & Inp,  int &Line, std::string PreviousSubset="", std::string AxisCategory="Axis",std::string AlternateIndex="");
    std::string strNode(int Level=0) const;
    std::vector<Axis> toVector() const;
    std::string subset() const {return _subset;}
    static std::string readSubset(ReadInput & Inp, std::string AxisCategory, int Line, std::string Default);

    AxisTree * factor(std::vector<std::string> Factor) const;
    AxisTree * complement(std::vector<std::string> Factor) const;
    void nodeCopy(const AxisTree *Node, bool View);
    void print() const;
};

#endif // AXISTREE_H
