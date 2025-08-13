// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "purgeIndex.h"
#include "index.h"
#include "operatorTree.h"
#include "printOutput.h"
#include "basisSub.h"

PurgeIndex::PurgeIndex(Index* Idx): idx(Idx){}

PurgeIndex& PurgeIndex::usingOperatorTree(const OperatorTree* Op){
    if(idx == Op->iIndex){
        opTrees.push_back({ Op, true });
    }else if(idx == Op->jIndex){
        opTrees.push_back({ Op, false });
    }else{
        ABORT("Index given is not one of OperatorTree given");
    }

    return *this;
}

void PurgeIndex::registerUsed(const OperatorTree* Op, bool IIndex){
    if(IIndex)
        used.insert(Op->iIndex);
    else
        used.insert(Op->jIndex);
    for(int i=0; i<Op->childSize(); i++){
        registerUsed(Op->child(i), IIndex);
    }
}

bool PurgeIndex::purge(Index* Idx){
    if(Idx->hasFloor()) return true;

    std::vector<int> removed;
    int total = Idx->basis()->size();
    for(int i=Idx->childSize()-1; i>=0; i--){
        if(used.find(Idx->child(i)) == used.end()){
            Idx->childErase(i);
            removed.push_back(i);
        }else{
            if(!purge(Idx->child(i))){
                Idx->childErase(i);
                removed.push_back(i);
            }
        }
    }

    if(removed.size() > 0
            and removed.size() < total
            and Idx->basis()){
        Idx->setBasis(Idx->basis()->remove(removed));
    }
    return removed.size() < total;
}

void PurgeIndex::prepare(){
    for(auto v: opTrees) registerUsed(v.first, v.second);
}

bool PurgeIndex::isUsed(const Index* Idx) const{
    return used.find(Idx) != used.end();
}

void PurgeIndex::run(){
    int sizeBefore = idx->sizeCompute();

    if(used.size() == 0) prepare();
    purge(idx);

    int sizeAfter = idx->sizeCompute();

    PrintOutput::DEVmessage("PurgeIndex: "+std::to_string(sizeBefore)+" -> "+std::to_string(sizeAfter));
}
