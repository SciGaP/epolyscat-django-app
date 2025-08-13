// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "overlapDVR.h"

#include "coefficients.h"
#include "coefficientsLocal.h"
#include "index.h"
#include "mpiWrapper.h"
#include "parallelOperator.h"
#include "operatorZD.h"
#include "parallelOperator.h"

using namespace std;

static void postDVR(OverlapDVR* Ovr){
    Ovr->postProcess();
    ParallelOperator::bcast(Ovr);
}

OverlapDVR::OverlapDVR(const Index *Idx)
    :OperatorTree("Ovr("+Idx->hierarchy()+")",Idx,Idx)
{
    // extract diagonal into Coefficients structure, force check
    Coefficients diagonal(iIndex);
    diagonal=Idx->localOverlap()->diagonal(true);
    _construct(&diagonal);
}
OverlapDVR::OverlapDVR(const OperatorTree *OvTree)
    :OperatorTree("Ovr("+OvTree->iIndex->hierarchy()+")",OvTree->iIndex,OvTree->jIndex)
{   
    if(not OvTree->isDiagonal())DEVABORT(OvTree->str()+"\nOverlapDVR only from diagonal Overlap");
    Coefficients diagonal(iIndex);
    diagonal=OvTree->diagonal(true);
    _construct(&diagonal);
}

OverlapDVR::OverlapDVR(Coefficients* Diagonal)
    :OperatorTree("Ovr("+Diagonal->idx()->hierarchy()+")",Diagonal->idx(),Diagonal->idx()){
    _construct(Diagonal);
}

void OverlapDVR::_construct(Coefficients* Diagonal){
    if(!Diagonal->parent())Diagonal->makeContinuous();
    if(Diagonal->idx()->hasFloor()){
        Eigen::MatrixXcd m(Diagonal->size(),Diagonal->size());
        for(int k=0;k<Diagonal->size();k++)m(k,k)=Diagonal->data()[k];
        oFloor=new OperatorZD(m,"ovr");
    }
    for(int k=0;k<Diagonal->childSize();k++)childAdd(new OverlapDVR(Diagonal->child(k)));
    if(!Diagonal->parent())postDVR(this);
}
