// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "indexProd.h"

using namespace std;
IndexProd::IndexProd(const Index *A, const Index *B)
{
    nodeCopy(A,false);
    for(size_t k=0;k<A->childSize();k++)
        childAdd(A->child(k)->deepCopy());

    Index * next=firstLeaf();
    while(next!=0){
        Index * node=next;
        next=next->nextLeaf();

        node->nodeCopy(B,false);
        for(size_t k=0;k<B->childSize();k++){
           node->childAdd(B->child(k)->deepCopy());
        }
    }
    if(isRoot())sizeCompute();
}
