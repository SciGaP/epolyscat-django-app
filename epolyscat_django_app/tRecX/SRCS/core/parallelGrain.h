// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARALLELGRAIN_H
#define PARALLELGRAIN_H

#include <vector>
#include "derivativeBlock.h"
#include <iostream>
#include <coefficients.h>

class Index;
class ParallelProcess;
/** \ingroup Parallelization */
/// blocks with equal left- and right hand Index's
class ParallelGrain
{
    friend class Parallel;
    friend class ParallelProcess;
    friend class DerivativeFlat;

    const ParallelProcess * leftProc; // owner process of left hand index, set after ParallelProcess::indexOwner
    std::vector<DerivativeBlock*> block;
public:
    ParallelGrain():leftProc(0){}
    inline const Index* leftIndex() const  {return block[0]->cInOut[1]->idx();}
    inline const Index* rightIndex() const {return block[0]->cInOut[0]->idx();}

    void setTargetProc(); // to be called after ParallelProcess::indexOwner has been set up
    double load() const {
        double l=0.;
        for(std::vector<DerivativeBlock*>::const_iterator b=block.begin();b!=block.end();b++)l+=(*b)->load();
        return l;
    }
    std::string str() const;
};

#endif // PARALLELGRAIN_H
