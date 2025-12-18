// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorFProduct.h"


OperatorFProduct::~OperatorFProduct(){
    if(not _view)for(auto f: _floors)delete f;
}

OperatorFProduct::OperatorFProduct(std::complex<double> Factor, std::vector<const OperatorFloor*> Floors, bool View)
    :OperatorFloor(Floors.front()->rows(),Floors.back()->cols(),"FProduct"),_factor(Factor),_floors(Floors),_view(View)
{
    if(_floors.size()<2)DEVABORT("need at least two floors fro OperatorFProduct");
    for(size_t i=0;i<_floors.size()-1;i++)
        _tmps.push_back(std::vector<std::complex<double>>(_floors[i]->cols()));
}

void OperatorFProduct::axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
                            const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const{
    _floors.back()->apply(Alfa*_factor,X,SizX,0.,_tmps.back().data(),_tmps.back().size());
    for(int i=_floors.size()-2;i>0;i--){
        _floors[i]->apply(1.,_tmps[i].data(),_tmps[i].size(),0.,_tmps[i-1].data(),_tmps[i-1].size());
    }
    _floors.front()->apply(1.,_tmps.front().data(),_tmps.front().size(),Beta,Y,SizY);
}

void OperatorFProduct::pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const{
    DEVABORT("not implemented");
}
