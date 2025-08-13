// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "indexExtract.h"

#include "mpiWrapper.h"
#include "printOutput.h"
#include "threads.h"
#include "indexConstraint.h"
#include "basisSub.h"
#include "overlapDVR.h"
#include "inverseDvr.h"

static std::vector<std::string> ExtractAxes;


//NOTE: leaf's must always be accepted
static bool selectAxes(const Index* Idx){
    return std::find(ExtractAxes.begin(),ExtractAxes.end(),Idx->axisName())!=ExtractAxes.end();
}

Index* IndexExtract::get(const Index *Idx, const std::vector<std::string> Ax, bool Complement, const IndexConstraint *Constraint){
    // count expected levels
    size_t levels=0;
    std::vector<std::string> hier(tools::splitString(Idx->hierarchy(),'.'));
    for(std::string ax: Ax)levels+=(std::find(hier.begin(),hier.end(),ax)==hier.end())==Complement;
    if(levels<Ax.size())return 0;
    return new IndexExtract(Idx,Ax,Complement,Constraint);
}

// first node that matches Select starting from Idx
static const Index* firstMatch(const Index *Idx, selectIndex Select){
    const Index* idx=Idx;
    // Note: leaf's cannot be omitted - always match
    while(idx!=0 and not (idx->isLeaf() and idx->axisName()=="NONE") and not Select(idx)){
        idx=idx->nodeNext(Idx);
    }
    return idx;
}

IndexExtract::IndexExtract(const Index *Idx, const std::vector<std::string> Ax, bool Complement, const IndexConstraint *Constraint, bool Entry)
{

    size_t levels=0;
    if(Entry){
        Idx->bottomExpandAll(); // make sure bottom leaf's are there

        // use axes or complement
        std::vector<std::string> hier(tools::splitString(Idx->hierarchy(),'.'));
        if(not Complement)
            ExtractAxes=Ax;
        else {
            // all axes NOT in Ax (avoid duplications)
            ExtractAxes.clear();
            for(std::string ax: hier)
                if(std::find(Ax.begin(),Ax.end(),ax)==Ax.end() and
                        std::find(ExtractAxes.begin(),ExtractAxes.end(),ax)==ExtractAxes.end())
                    ExtractAxes.push_back(ax);
        }
        // count expected levels
        for(std::string ax: ExtractAxes)levels+=std::count(hier.begin(),hier.end(),ax);
        if(levels<ExtractAxes.size())return;
    }

    const Index* idx=firstMatch(Idx,selectAxes);
    if(not idx)DEVABORT("algorithm error");
    setAxisName(idx->axisName());
    setBasis(BasisSub::superBas(idx->basis()));
    if(idx->depthInFloor()!=Index::npos)setFloor(0);

    // attach placeholder dummy leafs
    for(size_t k=0;k<basis()->size();k++)childAdd(new Index());
    std::vector<int>subset;
    // loop through first extract level within current subtree
    for(;idx!=0;idx=idx->nodeRight(Idx->parent()?Idx:0)){
        // basis must be (subsets of) the same basisAbstract
        if(not (*BasisSub::superBas(idx->basis())==*basis()))
            DEVABORT(Sstr+"bases on axis"+idx->axisName()+"differ within subtree"+Idx->index()+":"+Idx->strNode());

        // basis index as function of branch-index
        std::vector<int> basN=BasisSub::subset(idx->basis());

        // collect subtrees from all all branches of given idx
        for(const Index* bra=idx->descend();bra!=0;bra=bra->nodeRight(idx)){
            size_t nFun=basN[bra->nSibling()];
            // if not occpied yet, attach branch subtree
            if(not child(nFun)->basis()){
                subset.push_back(nFun);
                childReplace(nFun,new IndexExtract(bra,Ax,Complement,Constraint,false));
            }
        }
        if(subset.size()==basis()->size())break; // all filled, cut short
    }

    // remove placeholder branches
    // NOTE: handling of zero children is hacky - needs improvement
    if(subset!=BasisSub::subset(basis())){
        for(size_t k=childSize();k>0;k--)
            if(not child(k-1)->basis())childErase(k-1);
        setBasis(BasisAbstract::factory(BasisSub::strDefinition(basis(),subset)));
    }
    else
        setBasis(BasisSub::superBas(basis()));


    // post process
    if(Entry){
        if(Constraint!=0)Constraint->apply(this,{},{});
        cleanFloors();
        Idx->bottomUnexpandAll();
        sizeCompute();

        // this should go into Index member
        localOverlapAndInverse(0,0);
        IndexOverlap::set(this);
        Inverse::factory(this);
        bottomUnexpandAll();
    }

}

