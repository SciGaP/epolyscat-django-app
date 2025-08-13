// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "averageOverAxes.h"

#include "index.h"


void AverageOverAxes::setupIndex(const Index* root, Index* parent){
    if(not averageAxLevel[root->depth()]){
        for(int i=0; i<root->childSize(); i++){
            Index* self = new Index();
            self->nodeCopy(root->child(i), true); // View==true good idea?
            parent->childAdd(self);
            setupIndex(root->child(i), self);
        }
    }else{
        for(int i=0; i<root->childSize(); i++){
            setupIndex(root->child(i), parent);
        }
    }
}
    
void AverageOverAxes::averageRec(const Coefficients& src1, const Coefficients& src2, Coefficients& target){
    if(src1.idx() != src2.idx()) ABORT("Mismatch");

    if(src1.idx()->hasFloor()){
        if(not target.idx()->hasFloor()) ABORT("Internal");

        std::vector<std::complex<double>* > ptarget;
        target.pointerToC(ptarget);

        if(ptarget.size()!=1) ABORT("Internal");

        *(ptarget[0]) += src1.innerProductUnscaled(&src2);
        
    }else{
        if(averageAxLevel[src1.depth()]){
            for(int i=0; i<src1.childSize(); i++){
                averageRec(*src1.child(i), *src2.child(i), target);
            }
        }else{
            if(src1.childSize()!=target.childSize()) ABORT("Not tensor product structure");
            for(int i=0; i<src1.childSize(); i++){
                averageRec(*src1.child(i), *src2.child(i), *target.child(i));
            }
        }
    }
}

AverageOverAxes::AverageOverAxes(const Index* idx, std::vector<std::string> axNames): temp(idx){

    for(const Index* tmp=idx; tmp!=0; tmp=tmp->descend()){
        averageAxLevel.push_back(false);
        for(std::string ax: axNames){
            if(tmp->axisName().find(ax)!=std::string::npos){
                averageAxLevel.back()=true;
                break;
            }
        }

        if(tmp->hasFloor() and not averageAxLevel.back()) ABORT("Not supported!");
    }
    

    Index* Idx = new Index();
    Idx->nodeCopy(idx, true); // View==true good idea?
    setupIndex(idx, Idx);
    Idx->setFloor(Idx->firstLeaf()->depth());
    Idx->sizeCompute();

    this->Idx = Idx;
}

void AverageOverAxes::average(const Coefficients& src, Coefficients& target){
    target.setToZero();

    if(src.idx()->overlap()) src.idx()->overlap()->apply(1., src, 0., temp);
    else temp=src;

    averageRec(src, temp, target); 
}







