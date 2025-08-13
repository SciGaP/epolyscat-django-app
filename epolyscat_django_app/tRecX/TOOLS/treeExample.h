// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef TREEEXAMPLE_H
#define TREEEXAMPLE_H

#include "tree.h"

/// test class for tree
class TreeExample:public Tree<TreeExample>
{
public:
    int val;
    TreeExample():val(0){}
    TreeExample(unsigned int Depth, unsigned int Width, unsigned int LevelVal=1, const TreeExample* Parent=0);
    void valIndex(); // set value to index

    std::string strNode(int Precision=0) const {return tools::str(val);}

    bool nodeEquivalent(const TreeExample *Node) const{return val/100000==Node->val/100000;} //< this is not really a good example, as rearrangement would not be allowed
    bool nodeEmpty() const {return val<100000;}
    void nodeCopy(const TreeExample *Node,bool View){val=Node->val;}

    static void test();
};


#endif // TREEEXAMPLE_H
