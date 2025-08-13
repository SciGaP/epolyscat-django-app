// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coefficientsPermute.h"

#include "index.h"
#include "coefficientsFloor.h"
#include "printOutput.h"

using namespace std;

CoefficientsPermute::CoefficientsPermute(const Index *IOrig,std::vector<unsigned int> Perm)
{

    _cData=0;
    if(IOrig==0)return;

    // remove trivial permuations from tree floor
    while(Perm.size()>0 and Perm.back()==Perm.size()-1)Perm.pop_back();

    if(Perm.size()>0){
        Index* iOrig=IOrig->deepCopy();
        int levelOrig=iOrig->heightAboveFloor();
        iOrig->resetFloor(max(int(Perm.size())+1,levelOrig));

        // get permuted index
        Index * iNew=new Index();
        iOrig->permute(Perm,*iNew,false);
        iNew->resetFloor(max(int(Perm.size()),levelOrig));
        setIdx(iNew);
        delete iOrig;
    } else {
        setIdx(IOrig);
    }
    idx()->sizeCompute();

    if (idx()->hasFloor() or idx()->isLeaf()) {
        storageAssign(size(),0.);
    } else {
        // while above floor, continue descend
        for (unsigned int k=0;k<idx()->childSize();k++)
            childAdd(new CoefficientsPermute(idx()->child(k)));
    }

//    // storage (new default is to have all storage contigous (except in views)
    treeOrderStorage();

    // back-permuted view of this (i.e. view as the original)
    vector<unsigned int> back(Perm.size());
    for(int k=0;k<Perm.size();k++)back[Perm[k]]=k;
    permute(back,_viewPermAsOrig,true);
}

CoefficientsPermute & CoefficientsPermute::fromOrig(const Coefficients &COrig){
    copyView(true,const_cast<Coefficients*>(&COrig),&_viewPermAsOrig,0);
    return *this;
}

Coefficients & CoefficientsPermute::toOrig(Coefficients & COrig) const {
    copyView(false,&COrig,const_cast<Coefficients*>(&_viewPermAsOrig),0);
    return COrig;
}

void CoefficientsPermute::copyView(bool FromOrig, Coefficients *COrig, Coefficients *View, complex<double>*BData) const {
    if(COrig->isLeaf()){
        if(COrig->hasFloorData())BData=COrig->data();
        // data is contiguous at orginal floor, view may be descend below original floor
        if(View->isLeaf()){
            if(FromOrig)memcpy(View->data(),BData,View->idx()->sizeStored()*sizeof(*View->data()));
            else        memcpy(BData,View->data(),COrig->idx()->sizeStored()*sizeof(*COrig->data()));
        }
        else
            for (unsigned int k=0; k<View->childSize(); k++){
                copyView(FromOrig,COrig,View->child(k),BData);
                BData+=View->child(k)->idx()->sizeStored(); // advance pointer to original data
            }
    } else {
        for (unsigned int k=0; k<View->childSize(); k++){
            copyView(FromOrig,COrig->child(k),View->child(k),0);
        }
    }
}


