// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "bandedOverlap.h"
#include "operatorFloor.h"
#include "index.h"

BandedOverlap::BandedOverlap(Index *Idx, unsigned int SubD, unsigned int SuperD):OperatorTree("Overlap",Idx,Idx){
    if(Idx->hasFloor()){
        UseMatrix mat;
        Idx->overlap()->matrix(mat);
        // "mat.reband" requires "friend class IndexCoulX;" in useMatrix.h!
        UseMatrix *bandMat = new UseMatrix(mat.reband(SubD,SuperD,true));
        floor() = OperatorFloor::factory(std::vector<const UseMatrix*>(1,bandMat),"BandOvr"+Idx->hash()+Idx->hash());
    }
    else{
        for(int k=0;k<Idx->childSize();k++)
            childAdd(new BandedOverlap(Idx->child(k),SubD,SuperD));
    }
}
