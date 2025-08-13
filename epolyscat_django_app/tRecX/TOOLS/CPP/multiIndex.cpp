// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "multiIndex.h"

MultiIndex::MultiIndex(const std::vector<int> &M):m(M){}

void MultiIndex::first(std::vector<int> &i){i.assign(m.size(),0);}

/// increment multi-index, initial=empty vector, final=return false and empty vector
/// NOTE: rightmost indices run fastest
bool MultiIndex::next(std::vector<int> & i){
    if(i.size()==0){
        if(m.size()==0)return false;
        i.assign(m.size(),0);
        return true;
    }
    for(int d=m.size()-1;d>=0;d--){
        if(i[d]<m[d]-1){
            i[d]++;
            for(unsigned int l=d+1;l<m.size();l++)i[l]=0;
            return true;
        }
    }
    // if it cannot be incremented, reset and return false
    i.clear();
    return false;
}
bool MultiIndex::nextCol(std::vector<int> & i){
    if(i.size()==0){
        if(m.size()==0)return false;
        i.assign(m.size(),0);
        return true;
    }
    for(size_t d=0;d<m.size();d++){
        if(i[d]<m[d]-1){
            i[d]++;
            for(size_t l=0;l<d;l++)i[l]=0;
            return true;
        }
    }
    // if it cannot be incremented, reset and return false
    i.clear();
    return false;
}
