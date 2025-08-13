// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORIDENTITY_H
#define OPERATORIDENTITY_H

#include "operatorTree.h"
#include "operatorDefinition.h"

#include "index.h"
#include "parallel.h"
#include "operatorFloor.h"

/// tree-structured identity operator
///
/// when tree-structure is not needed use OperatorIdentity for maximal efficiency
class OperatorIdentityTree : public OperatorTree
{
    std::complex<double> _factor;
    OperatorIdentityTree(const Index* Idx,std::complex<double> Factor,bool Top)
        :OperatorTree("Identity",Idx,Idx){
        if(Idx->hasFloor()){
            oFloor=Parallel::operatorFloor(idx(),idx(),[&](){ return OperatorFloor::factory("<Id>",iIndex,jIndex,Factor);});
        }
        else {
            for(size_t k=0;k<Idx->childSize();k++)
                childAdd(new OperatorIdentityTree(Idx->child(k),Factor,false));
        }
        if(Top)postProcess();
    }
public:
    OperatorIdentityTree(const Index* Idx,std::complex<double> Factor=1.)
        :OperatorIdentityTree(Idx,Factor,true){}
};

/// direct identity operator
///
/// when tree-structure is needed use OperatorIdentityTree (with minor performance penalty)
class OperatorIdentity : public OperatorTree
{
    std::complex<double> _factor;
public:
    OperatorIdentity(const Index* Idx,std::complex<double> Factor=1.):OperatorTree("Identity",Idx,Idx),_factor(Factor){}
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
    {Y.axpy(A*_factor,Vec,B);}
};

#endif // OPERATORIDENTITY_H
