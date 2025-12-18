// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "diagnose.h"

#include "operatorAbstract.h"
#include "useMatrix.h"
#include "index.h"
#include "str.h"
#include "printOutput.h"

using namespace std;

Diagnose::Diagnose(){}

void Diagnose::show(const OperatorAbstract *Op,unsigned int BlockLevel) const {

    int iHeight=max(min(BlockLevel,Op->iIndex->height())-1,(unsigned int)(0));
    int jHeight=max(min(BlockLevel,Op->jIndex->height())-1,(unsigned int)(0));

    // this should rather be: Op->iIndex->Tree::size(iHeight);
    int matI=Op->iIndex->descend(iHeight)->levelSize();
    int matJ=Op->jIndex->descend(jHeight)->levelSize();

    if((matI+matJ)>200){
        PrintOutput::warning(Str("too large to show, choose smaller block level, is")+BlockLevel);
        return;
    }

    UseMatrix mat;
    Op->matrix(mat);

    UseMatrix showMat;
    if(matI!=Op->iIndex->sizeStored() or matJ!=Op->jIndex->sizeStored()){
        showMat=UseMatrix::Zero(matI,matJ);
        const Index * jdx=Op->jIndex->descend(jHeight);
        for(int j=0;j<matJ;j++,jdx=jdx->nodeNext()){
            const Index * idx=Op->iIndex->descend(iHeight);
            for(int i=0;i<matJ;i++,idx=idx->nodeNext()){
                showMat(i,j)=mat.block(idx->posIndex(Op->iIndex),jdx->posIndex(Op->jIndex),idx->sizeStored(),jdx->sizeStored()).maxAbsVal();
            }
        }
    }
    else
        Op->matrix(showMat);

    showMat.print(Str(Op->name)+"block ("+Op->iIndex->index()+"|"+Op->jIndex->index()+") to level "+BlockLevel+"hermitian:"+mat.isHermitian(),0);
}

