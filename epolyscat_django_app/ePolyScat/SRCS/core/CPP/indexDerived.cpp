// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "indexDerived.h"

using namespace std;

IndexDerived::IndexDerived(){}

string IndexDerived::strNode() const{
    string s=Index::strNode();
    if(fromIndex.size()!=0){
        s+=" from: "+tools::str(fromIndex[0]->index());
        for(unsigned int k=1;k<fromIndex.size();k++)s+=" |"+tools::str(fromIndex[k]->index());
    }
    return s;
}
