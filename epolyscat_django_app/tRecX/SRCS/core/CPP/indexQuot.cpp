// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "indexQuot.h"

#include <string>

using namespace std;

IndexQuot::IndexQuot(const Index *Full, const Index *Fac)
{
    if(Full->isRoot()){
        // check
        for(int k=0;Fac->descend(k)!=0;k++)
            if(Fac->descend(k)->axisLevel(Full)==npos){
                cout<<"Fac\n"<<Fac->str()<<endl;
                cout<<"Full\n"<<Full->str()<<endl;
                ABORT("not all levels of Fac in Full");
            }
    }

    while(Full->axisLevel(Fac)!=npos and not Full->isLeaf()){
        Full=Full->descend();
        if(Full==0)return;
    }
    nodeCopy(Full,false);
    for (int k=0;k<Full->childSize();k++)
        childAdd(new IndexQuot(Full->child(k),Fac));

    if(Full->isRoot())sizeCompute();
}
