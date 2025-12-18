// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coefficientsViewDeep.h"

#include "index.h"
#include "coefficientsFloor.h"

CoefficientsViewDeep::CoefficientsViewDeep(const Index *Idx, int Depth, bool paraSub){
    _view.reset(Idx);

    _view.idx()->sizeCompute(); // make sure we have correct numbering

    // expand all leaf's to full depth
    Coefficients* f=_view.firstLeaf();
    for(;f!=0;f=f->nodeNext()){
        f->_cData=0;
        for(size_t k=0;k<f->idx()->childSize();k++)
            f->childAdd(new Coefficients(f->idx()->child(k)));
    }

    // truncate to depth
    f=_view.descend(std::min(Depth,int(Idx->height())));
    while(f!=0){
        while(f->childSize()>0)f->childErase(0);
        f->_cData=0;
        f=f->nodeNext();
    }
    disown(&_view); // abandon ownership
}

void CoefficientsViewDeep::disown(Coefficients*C){
    C->makeView();
    C->nodeStorageClear();
    for(size_t k=0;k<C->childSize();k++)disown(C->child(k));
}

Coefficients* CoefficientsViewDeep::view(Coefficients* C){
    if(_view.idx()!=C->idx())ABORT("view does not match coefficients");
    view(&_view,C);
    return &_view;
}


void CoefficientsViewDeep::view(Coefficients* View, Coefficients* C){
    if(C->isLeaf())extend(View,C->floorData());
    for(size_t k=0;k<C->childSize();k++)
        view(View->child(k),C->child(k));
}

void CoefficientsViewDeep::extend(Coefficients* View, std::complex<double> * Data){
    if(View->isLeaf()){
        View->setFloorData(Data);
    }
    for(size_t k=0;k<View->childSize();k++){
        extend(View->child(k),Data);
        Data+=View->child(k)->idx()->sizeStored();
    }
}
