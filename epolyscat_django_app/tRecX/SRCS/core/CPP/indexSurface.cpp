// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "indexSurface.h"
#include "basisGrid.h"
#include "basisSub.h"
#include "basisAbstract.h"
#include "basisIntegrable.h"

using namespace std;

IndexSurface::IndexSurface(const Index *From, std::vector<double> Radius, unsigned int NSurf)
{
    if(From->isLeaf())ABORT(From->root()->str()+"\n\nNo unbounded axis found in Index");

    Index::nodeCopy(From,false); // true copy, not view
    fromIndex.push_back(From);

    // second occurance of name indicates multi-piece axis
    const Index* fIndex=From->descend()->axisIndex(From->axisName());

    if(fIndex==0 or NSurf>0)
    {   // not desired surface - descend
        if(0!=fIndex)NSurf--;
        for(unsigned int k=0;k<From->childSize();k++)childAdd(new IndexSurface(From->child(k),Radius,NSurf));
    }

    else
    {
        // indicate surface level
        setAxisName("surf"+From->axisName());
        fromIndex.push_back(From);
        vector<int> elemN;
        for(unsigned int k=0;k<Radius.size();k++)
        {
            // locate element containing Radius
            const Index* elem=From->descend();
            unsigned int toFunc=fIndex->depth()-elem->depth();

            for(;elem!=0;elem=elem->rightSibling()){

                Index* elemFunc=elem->descend(toFunc); // descend to function level
                double lB=elemFunc->basis()->integrable()->lowBound(),uB=elemFunc->basis()->integrable()->upBound();
                // locate inside element; if on boundary, choose element nearer to 0
                if(     (lB>=0 and lB< Radius[k] and Radius[k]<=uB) or
                        (lB< 0 and lB<=Radius[k] and Radius[k]< uB)    )
                {
                    // attach copy of tree that contains surface
                    for(unsigned int l=childSize();l>0;l--)childErase(l-1);
                    elemN.push_back(elem->nSibling());
                    childAdd(new Index(*From->child(elemN.back())));

                    // insert v/d at matching level
                    for(Index *vd=childBack()->descend(toFunc);vd!=0;vd=vd->nodeRight(),elemFunc=elemFunc->nodeRight()){
                        // erase function level
                        for(unsigned int l=vd->childSize();l>0;l--)vd->childErase(l-1);
                        // insert v/d level
                        // NOTE: we need all possible levels below to be equivalent
                        for(unsigned int l=0;l<2;l++){
                            vd->childAdd(new Index(*elemFunc->child(0)));
                                vd->child(l)->setBasis(BasisGrid::factory(std::vector<double>(1,Radius[k])));
                            vd->child(l)->leafAdd();
                        }

                        vd->setBasis(BasisAbstract::factory("Vector:2"));
                        vd->setAxisName("v/d");
                    }
                    break;
                }
            }
            if(elem==0)ABORT("surface outside all elements: "+tools::str(Radius[k]));
        }
        // basis on surf-level is a grid of element numbers
        setBasis(BasisAbstract::factory(BasisSub::strDefinition(From->basis(),elemN)));
    }

    if(From->depth()==0){
        // set floor to below v/d level (or to From floor, if deeper)
        resetFloor(max(From->firstFloor()->depth(),axisIndex("v/d")->depth()+1));
    }
}
