// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARALLELCROSS_H
#define PARALLELCROSS_H

#include <vector>
#include "derivativeBlock.h"
#include <iostream>
#include <string>
#include "mpiWrapper.h"

class Index;
class ParallelProcess;
/** \ingroup Parallelization */
/// rows and columns to be owned by one process
class ParallelCross
{
    friend class Parallel;
    friend class ParallelProcess;
    friend class DerivativeFlat;
    std::vector<DerivativeBlock*> colBlock,rowBlock;
public:
    ParallelCross(){}
    std::string str() const;
    const Index* index() const;
    const Index* nonOwnerIndex() const;
    double load(bool Sync=true) const {
        double l=0.;
        for(std::vector<DerivativeBlock*>::const_iterator b=rowBlock.begin();b!=rowBlock.end();b++)l+=(*b)->load();
        for(std::vector<DerivativeBlock*>::const_iterator b=colBlock.begin();b!=colBlock.end();b++)l+=(*b)->load();
        if(Sync)MPIwrapper::Bcast(&l,1,MPIwrapper::master());
        return l;
    }
};

#endif // PARALLELSTRIPE_H
