// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "parallelBuffer.h"
#include "coefficients.h"
#include <numeric>
#include "str.h"
#include "timer.h"
void ParallelBuffer::add(const Coefficients *C){
    if(C->isContiguous())
        add(C->anyData(),C->size());
    else
        for(int k=0;k<C->childSize();k++)add(C->child(k));

}

int ParallelBuffer::extract(int Pos,Coefficients *C){
    int pos(Pos);
    if(C->isContiguous())
        pos=extract(pos,C->anyData(),C->size());
    else
        for(int k=0;k<C->childSize();k++)pos=extract(pos,C->child(k));
    return pos;
}

void ParallelBuffer::allGather(){
    std::vector<int> gsiz(MPIwrapper::Size(),0);
    gsiz[MPIwrapper::Rank()]=val.size();
    MPIwrapper::AllreduceSUM(gsiz.data(),gsiz.size());// time-consuming -not a large amount when multiple nodes
    std::vector<std::complex<double> > gval(std::accumulate(gsiz.begin(),gsiz.end(),0));
    if(gval.size()==0)return; //empty
    MPIwrapper::AllGatherV(val.data(),val.size(),gval.data(),gsiz.data());
    std::swap(val,gval);
}
