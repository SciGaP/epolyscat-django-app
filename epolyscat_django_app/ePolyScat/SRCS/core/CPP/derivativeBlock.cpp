// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "derivativeBlock.h"

#include "coefficients.h"
#include "useMatrix.h"
#include "operatorTree.h"
#include "operatorFloor.h"
#include "str.h"

using namespace std;


DerivativeBlock::DerivativeBlock(std::vector<unsigned int> BlockSort,
                                 const OperatorTree* OLeaf, Coefficients * ICoef, Coefficients*JCoef, double ONorm, double *XNorm)
    :oLeaf(OLeaf),xNorm(XNorm),blockSort(BlockSort){

    eps=ONorm;
    cInOut.push_back(JCoef);
    cInOut.push_back(ICoef);
#ifdef _OPENMP
    _alt=0;
#endif

    // distribute the send/receive info
    // this should be replaced by ParallelOperator
    vector<complex<double> > buf;
}

double DerivativeBlock::load() const{return oLeaf->floor()->applicationCost(false);}

bool DerivativeBlock::lessEqual(const DerivativeBlock &A, const DerivativeBlock &B){
    for(unsigned int k=0;k<A.blockSort.size();k++)
        if(A.blockSort[k]!=B.blockSort[k])
            return A.blockSort[k]<B.blockSort[k];
    return &A != &B;
}

bool DerivativeBlock::less(const DerivativeBlock &A, const DerivativeBlock &B){
    for(unsigned int k=0;k<A.blockSort.size();k++)
        return A.blockSort[k]<B.blockSort[k];
    return true;
}
