// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "indexPermute.h"

IndexPermute::IndexPermute(const Index *Idx, std::string Hierarchy)
{
    // get permutation from Idx.hierarchy to desired
    std::vector<std::string> axP=tools::splitString(Hierarchy,'.');
    std::vector<std::string> axI=tools::splitString(Idx->hierarchy(),'.');
    std::vector<unsigned int> perm(axI.size());
    for(size_t k=0;k<axP.size();k++){
        int cntP=std::count(axP.begin(),axP.begin()+k,axP[k]);
        for(size_t l=0;l<axI.size();l++)
            if(axP[k]==axI[l] and cntP==std::count(axI.begin(),axI.begin()+l,axI[l]))perm[l]=k;
    }
    // expand Idx to leaf level (if required)
    Idx->bottomExpandAll();
    Idx->permute(perm,*this,false);
    size_t floor=Idx->heightAboveFloor();
    sizeCompute();

    // create map
    resetFloor(heightAboveBottom()+1);
    const_cast<Index*>(Idx)->resetFloor(Idx->heightAboveBottom()+1);
    Coefficients c(Idx),p(this);
    std::unique_ptr<Coefficients> pViewC(c.permute(perm,true));
    for(size_t k=0;k<c.size();k++)c.data()[k]=k;
    p=*pViewC;
    for(size_t k=0;k<p.size();k++)_permC.push_back(size_t(p.data()[k].real()));

    // restore expansion state
    resetFloor(floor);
    bottomExpandAll();
    const_cast<Index*>(Idx)->resetFloor(floor);
    Idx->bottomUnexpandAll();

}
