// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorInverseCoulX.h"
#include "index.h"
#include "operatorFloor.h"

OperatorInverseCoulX::OperatorInverseCoulX(const Index *Idx, unsigned int SubD, unsigned int SuperD, bool BandOvr):OperatorTree("Inv(CoulXOvr)",Idx,Idx){
    if(Idx->hasFloor()){
        oFloor = OperatorFloor::factoryInverse(Idx, SubD, SuperD, BandOvr);
    }
    else{
        for(int k=0;k<Idx->childSize();k++){
            childAdd(new OperatorInverseCoulX(Idx->child(k),SubD,SuperD,BandOvr));
        }
    }
}
