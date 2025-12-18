// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "indexOverlap.h"

#include "basisAbstract.h"
#include "index.h"
#include "operatorTree.h"
#include "operatorDefinition.h"
#include "printOutput.h"
#include <typeinfo>
#include "overlapDVR.h"
#include "inverseFloors.h"
#include "mpiWrapper.h"
#include "parallelOperator.h"
#include "operatorIdentity.h"

struct OverlapAndInverse{
    std::vector<OperatorTree*> trees;
    Inverse* inverse;
};
static std::map<const Index*,OverlapAndInverse> _topOverlapList;
static OverlapAndInverse & getFromTopOverlap(const Index* Idx, std::vector<unsigned int> &idx, const Index *&iRoot){

    // descend from root towards node until top operator found
    idx=Idx->index();
    std::map<const Index*,OverlapAndInverse>::iterator p;

    iRoot=Idx->root();
    for(;0!=iRoot
        and _topOverlapList.end()==(p=_topOverlapList.find(iRoot))
        and iRoot->depth()<idx.size();
        iRoot=iRoot->child(idx[iRoot->depth()]));
    if(p==_topOverlapList.end())DEVABORT(Str("not found in _topOverlapList")+Idx->strNode(0)+idx);
    return p->second;
}

const OperatorTree* IndexOverlap::getOverlap(const Index* Idx,std::vector<unsigned int> &idx, const Index *&iRoot, int Kind){
    return getFromTopOverlap(Idx,idx,iRoot).trees[Kind];
}

OperatorTree * IndexOverlap::getTree(const Index* Idx, int Kind) {

    auto top=_topOverlapList.find(Idx);
    if(top!=_topOverlapList.end())return top->second.trees[Kind];

    std::string mess=Sstr+"overlap list failed at node ["+Idx+"] ( level"+Idx->depth()+") in hierachy"+Idx->root()->hierarchy()+
            "\npossible reasons:\n 1. attempted use in permuted index or other index w/o overlap\n 2. not implemented";

    const Index* idx=Idx;
    while(idx and (top=_topOverlapList.find(idx))==_topOverlapList.end())idx=idx->parent();
    if(!idx)DEVABORT(Sstr+"no entry found for \n"+Idx->str()+"at"+Idx+"\n"+mess);
    if(Idx==idx)return top->second.trees[Kind];

    // descend to present node
    OperatorTree* opNode=top->second.trees[Kind];
    // CAUTION: may be very slow
    for(;opNode;opNode=opNode->nodeNext()){
        if(opNode->iIndex==Idx and opNode->jIndex==Idx)break;
    }
    if(!opNode)
        DEVABORT(Sstr+"indices do not match for Kind="+Kind+"\nthis\n"+Idx+Idx->str()+"\nfound\n"+top->second.trees[Kind]->str()+"\n"+mess);

    return opNode;
}

// NOTE: we need overlap operators to be available everywhere
OperatorTree* bcastDiagonal(OperatorTree* Op){
    ParallelOperator::bcast(Op);
    return Op;
}

bool IndexOrthonormal(const Index* Idx){
    if(not Idx->basis()->isOrthonormal())return false;
    for(size_t k=0;k<Idx->childSize();k++)
        if(not IndexOrthonormal(Idx->child(k)))return false;
    return true;
}

void IndexOverlap::localOverlapAndInverse(const Index* Idx, OperatorTree* Ovr, OperatorTree* Inv){

    if(Idx->parent()==0 and (Ovr!=0 or Inv!=0))DEVABORT("must enter from root with Ovr=Inv=0");


    if(_topOverlapList.count(Idx->parent()))
        return; // already set up from higher up - nothing to be done

    // names containing "&" are of hybrid axes
    std::vector<std::string> parts=tools::splitString(Idx->axisName(),'&');
    if(parts.size()>1){
        // all blocks below are non-hybrid
        OperatorTree* ov=new OperatorTree("Overlap",Idx,Idx);
        for(size_t k=0;k<Idx->childSize();k++){
            if(parts[k]!=Idx->child(k)->axisName())
                PrintOutput::DEVwarning(Sstr+"incorrect construction of hybrid axis: "+Idx->axisName()+k+Idx->child(k)->axisName()
                                        +"\n"+Idx->root()->str());
            localOverlapAndInverse(Idx->child(k),0,0);
            ov->childAdd(Idx->child(k)->localOverlap());
        }
        if(not _topOverlapList.count(Idx->parent()))_topOverlapList[Idx].trees.push_back(ov);
    }

    else if(IndexOrthonormal(Idx) or Idx->hierarchy().find("Ndim.Ndim")==0){
        if(not _topOverlapList.count(Idx->parent())){
            _topOverlapList[Idx].trees.push_back(bcastDiagonal(new OperatorTree("Overlap","[[Id]]",Idx,Idx)));
        }
    }
    else if(Ovr==0){
        if(not _topOverlapList.count(Idx->parent())){
            OperatorTree* ov=new OperatorTree("Overlap","<<Overlap>>",Idx,Idx);
            _topOverlapList[Idx].trees.push_back(bcastDiagonal(ov));
        }
    }
    else {
        if(not _topOverlapList.count(Idx->parent()))_topOverlapList[Idx].trees={Ovr,Inv};
    }

    // get all diagonal floor blocks inverted
    if(_topOverlapList[Idx].trees.size()==1){
        _topOverlapList[Idx].trees.push_back(new InverseFloors(_topOverlapList[Idx].trees[0]));
    }

}
static std::map<const Index*,const OperatorAbstract*> _overlapList;
const OperatorAbstract* IndexOverlap::get(const Index *Idx){
    if(not _overlapList.count(Idx))return 0;
    return _overlapList[Idx];
}

void setAllBlocks(const OperatorAbstract* Ovr){
    if(Ovr->iIndex==Ovr->jIndex){
        _overlapList[Ovr->iIndex]=Ovr;
        const OperatorTree* ovTree=dynamic_cast<const OperatorTree*>(Ovr);
        if(ovTree){
            if(ovTree->isBlockDiagonal() or ovTree->iIndex->axisSubset()=="Subspace&Complement")
                for(size_t k=0;k<ovTree->childSize();k++)setAllBlocks(ovTree->child(k));
        }
    }
}

// set overlap pointer: if not set, seek up in hiearchy, or else construct new
void IndexOverlap::set(const Index *Idx, const OperatorAbstract* Ov){
    if(Ov){
        const OperatorTree* ovTree=dynamic_cast<const OperatorTree*>(Ov);
        if(ovTree->isDiagonal()){
            PrintOutput::DEVwarning("possible memory leak when recomputing Overlap "+Idx->strNode());
            Ov=new OverlapDVR(ovTree);
        }
        if(Idx->overlap()){
            // same already set
            if(Idx->overlap()==Ov)return;
            const OperatorTree* ov1=dynamic_cast<const OperatorTree*>(Ov);
            const OperatorTree* ov2=dynamic_cast<const OperatorTree*>(Idx->overlap());
            if(ov1 and ov2){
                if(Ov!=Idx->overlap())
                    PrintOutput::DEVwarning("re-compute of overlap - likely memory leak for "+Idx->strNode()+Idx->hierarchy());
                return;
            }
            PrintOutput::DEVwarning("must not reset overlap, Idx="+Idx->hierarchy());
        }
        setAllBlocks(Ov);
    }
    else {
        if(Idx->overlap())return;
        if(get(Idx))return;
        localOverlapAndInverse(Idx,0,0);

        if(!Idx->overlap() or !get(Idx)){
            // no old found - construct new
            if(IndexOrthonormal(Idx)){
                setAllBlocks(new OperatorIdentity(Idx));
            } else {
                OperatorDefinition def("<<Overlap>>",Idx->hierarchy());
                OperatorTree *ov0=new OperatorTree("Ovr("+Idx->hierarchy()+")","<<Overlap>>",Idx,Idx);
                ParallelOperator::bcast(ov0);
                if(ov0->isDiagonal()){
                    setAllBlocks(new OverlapDVR(ov0));
                    delete ov0;
                }
                else {
                    setAllBlocks(ov0);
                }
            }
        }
    }
}
static std::map<const Index*,const Inverse*> _inverseOverlapList;
void IndexOverlap::setInverseOverlap(const Index* Idx,const Inverse *Inv){
    if(Inv!=0)_inverseOverlapList[Idx]=Inv;
}
const Inverse * IndexOverlap::inverseOverlap(const Index* Idx) {
    if(_inverseOverlapList.count(Idx)==0)return 0;
    return _inverseOverlapList[Idx];
}
