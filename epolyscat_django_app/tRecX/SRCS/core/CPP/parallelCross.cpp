// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "parallelCross.h"
#include "coefficients.h"
#include "index.h"
#include "parallel.h"

const Index* ParallelCross::index() const
{
    if(colBlock.size()==0 and rowBlock.size()==0)DEVABORT("emtpy cross");
    return colBlock.size()?colBlock[0]->cInOut[0]->idx():rowBlock[0]->cInOut[1]->idx();
}
const Index* ParallelCross::nonOwnerIndex() const
{
    if(colBlock.size()==0 and rowBlock.size()==0)DEVABORT("emtpy cross");
    return colBlock.size()?colBlock[0]->cInOut[1]->idx():rowBlock[0]->cInOut[0]->idx();
}

// show row/column to which the cross belongs
std::string ParallelCross::str() const {
    std::string s;
    const Index* idx(0);
    if(colBlock.size()>0)     idx=colBlock[0]->cInOut[0]->idx();
    else if(rowBlock.size()>0)idx=rowBlock[0]->cInOut[1]->idx();
    if(colBlock.size()>0)s+="col: "+colBlock[0]->oLeaf->root()->name+colBlock[0]->oLeaf->strNode();
    if(rowBlock.size()>0)s+="row: "+rowBlock[0]->oLeaf->root()->name+rowBlock[0]->oLeaf->strNode();
    s+=tools::str(idx->posIndex(),7)+"["+idx->hash()+"]";
    s+=":"+tools::str(int(rowBlock.size()))+"/"+tools::str(int(colBlock.size()));
    s+=" load "+tools::str(load(),2);
    s+=" owner: "+index()->strNode()+" other: "+nonOwnerIndex()->strNode();
//    s+="\n"+Parallel::strDist();

    return s;
}



