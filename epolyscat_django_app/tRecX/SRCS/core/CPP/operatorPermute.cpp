// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorPermute.h"


#include "index.h"
#include "coefficients.h"
#include "basisVector.h"


OperatorPermute::OperatorPermute(std::string Name, size_t JLevel, const Index *JIndex, bool MapFromIndex)
    :OperatorTree(Name,0,0)
{
    if(JLevel<1)DEVABORT("cannot permute around 0-level");
    if(JLevel>JIndex->heightAboveBottom())
        DEVABORT(Sstr+"cannot permute at level "+JLevel+"- at or below bottom level"+JIndex->heightAboveBottom());

    // determine the permuted subindex sizes
    std::vector<int> iSiz(JIndex->descend(JLevel)->levelSize(),0);
    for(const Index* jT=JIndex->descend(JLevel);jT!=0;jT=jT->nodeRight(JIndex)){
        for(size_t i=0;i<iSiz.size();i++)
            if(jT->size()>i)iSiz[i]++;
    }

    // construct the permuted Index (should be moved into routine)
    Index *ix=new Index();
    std::vector<size_t>iPos(1,0);
    for(size_t n=0;n<iSiz.size();n++)
        if(iSiz[n]){
            iPos.push_back(iPos.back()+iSiz[n]);
            ix->childAdd(new Index({new BasisVector(iSiz[n])},{"Vec"}));
        }
    iPos.pop_back();
    ix->setBasis(new BasisVector(ix->childSize()));
    ix->setAxisName("Vec");
    ix->resetFloor(1);
    ix->sizeCompute();
    iIndex=ix;
    jIndex=JIndex;

    if(iIndex->size()!=JIndex->size())
        DEVABORT("perm"+iIndex->str(1)+"\n"+JIndex->str(1)
                 +"\npermuted index size does not match original");

    // get the permutation map
    _perm.resize(JIndex->size());
    size_t jPos=0;
    for(const Index* jT=JIndex->descend(JLevel);jT!=0;jT=jT->nodeRight(JIndex))
        for(size_t p=0;p<jT->size();iPos[p]++,jPos++,p++){
            if(MapFromIndex)_perm[jPos]=iPos[p];
            else            _perm[iPos[p]]=jPos;
        }


    if(not MapFromIndex)std::swap(iIndex,jIndex);
}

void OperatorPermute::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    if(iIndex!=Y.idx())DEVABORT("lhs index does not match");
    if(jIndex!=Vec.idx())DEVABORT("rhs index does not match");
    if(not Vec.orderedData() or not Y.orderedData())DEVABORT("need contiguous storage on permute level");
    for(size_t k=0;k<_perm.size();k++)Y.orderedData()[_perm[k]]=A*Vec.orderedData()[k]+B*Y.orderedData()[_perm[k]];
}


std::string OperatorPermute::strNode(int Digits) const{
    if(not isLeaf())return "";
    return tools::str(_perm,3,"");
}
