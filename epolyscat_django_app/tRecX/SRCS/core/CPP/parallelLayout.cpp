// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "parallelLayout.h"

#include <algorithm>

#include "index.h"
#include "parallelCross.h"
#include "parallel.h"
#include "readInput.h"
#include "printOutput.h"
#include "mpiWrapper.h"

static std::map<std::vector<std::string>,std::string> _listSorting;

/// set sorting for Coors (unless already one set)
static void setSort(std::string Coors){
    std::vector<std::string> vCoor=tools::splitString(Coors,'.');
    std::sort(vCoor.begin(),vCoor.end());
    _listSorting[vCoor]=Coors;
}

/// default sortings
static void setup(){
    setSort("X1.X2");
    setSort("Phi.Eta.Rn");
    setSort("Phi1.Phi2.Eta1.Eta2.Rn1.Rn2");
    setSort("Orbital&Phi.Phi.Eta.Rn");
    setSort("Subspace.SubSubspace");
}

void ParallelLayout::read(ReadInput &Inp){
    setup();
    std::string val;
    Inp.read("Parallel","sort",val,"Phi.Eta.Rn","coordinate sorting for parallelization, first runs fastest");
    setSort(val);
}

ParallelLayout::ParallelLayout(const Index* Idx, std::string Sorting)
{
    // get hierarchy and sorted hierarchy
    _hier=Idx->hierarchy();
    std::vector<std::string> vHier=tools::splitString(_hier,'.');
    // remove certain special axis names and duplicates
    for(int k=vHier.size()-1;k>=0;k--){
        if(        vHier[k].find("spec")==0
                   or vHier[k].find("Subspace")==0
                   or std::find(vHier.begin(),vHier.begin()+k,vHier[k])!=vHier.begin()+k
                   )
            vHier.erase(vHier.begin()+k);
    }
    std::string hierSort;
    for(std::string ax: vHier)hierSort+="."+ax;
    hierSort.erase(hierSort.begin());     // remove leading "."
    std::sort(vHier.begin(),vHier.end()); // bring to standard (lexicographical) sorting
    if(vHier.size()==1){
        _permute.assign(1,0);
        return; // no sorting needed
    }

    // seek predefined sorting
    if(Sorting==""){
        if(_listSorting.size()==0)setup();
        if(_listSorting.count(vHier))Sorting=_listSorting[vHier];
    }

    if(Sorting!=""){
        if(MPIwrapper::Size()>1)PrintOutput::message("using ParallelLayout "+Sorting);
        std::vector<std::string> vSort=tools::splitString(Sorting,'.');
        if(vHier.size()!=vSort.size())ABORT(Sstr+"index hierarchy does not match Sorting: "+vHier+" --- "+Sorting);
        for(size_t k=0;k<vSort.size();k++){
            if(std::find(vHier.begin(),vHier.end(),vSort[k])==vHier.end())
                ABORT(Sstr+"index hierarchy does not match Sorting: "+vHier+" vs "+Sorting);
            _permute.push_back(std::find(vHier.begin(),vHier.end(),vSort[k])-vHier.begin());
        }
    }
    else {
        if(MPIwrapper::Size()>1 and vHier.size()>1)
            PrintOutput::DEVwarning(Sstr+"no sorting defined for "+vHier
                                    +"\navailable\n"+tools::listMapKeys(_listSorting,"\n")
                                    +"\nusing hierarchy sorting, change through Parallel:sort");
        for(size_t k=0;k<Idx->height();k++)_permute.push_back(k);
    }
    if(Idx->root()->isHybrid()){
        PrintOutput::DEVwarning("no sorting for hybrid basis");
        _permute.clear();
    }
    else if(_permute.size()<2)
        PrintOutput::DEVwarning("permutation of size 1 - problem with algorithm, fixed by hack");
}

static std::vector<int> thisPerm; // locally communicate permutation into smallerCross
bool ParallelLayout::smallerCross(const ParallelCross * A, const ParallelCross * B){
    if(A->index()->index().size()!=B->index()->index().size()){
        // map between different size indices: let the larger index run faster
        return A->index()->root()->size()>B->index()->root()->size();
    }
    for(int k=A->index()->index().size();k>0;k--){
        size_t pk=k;
        if(A->index()->index().size()<=thisPerm.size()){
            pk=thisPerm.at(k-1);
            if(pk>=A->index()->index().size()){
                PrintOutput::DEVwarning(Sstr+"Perm outside size"+k+pk+thisPerm+thisPerm.size()+A->index()->strNode(),10);
                pk=std::min(A->index()->index().size(),B->index()->index().size())-1;
            }
        }
        if(A->index()->index()[pk]<=B->index()->index()[pk])return true;
        if(A->index()->index()[pk]> B->index()->index()[pk])return false;
    }
    return A->index()->size()>B->index()->size();
}

static bool smallerIndex(const Index * A, const Index * B){
    if(A->index().size()!=B->index().size()){
        // map between differen size indices: let the larger index run faster
        return A->root()->size()>B->root()->size();
    }
    for(int k=A->index().size();k>0;k--){
        int pk=thisPerm[k-1];
        if(pk>=A->index().size())DEVABORT(Sstr+"Perm outside size"+k+pk+thisPerm+thisPerm.size()+"\n"+A->strNode());
        if(A->index()[pk]<=B->index()[pk])return true;
        if(A->index()[pk]> B->index()[pk])return false;
    }
    return A>B;
}

int ParallelLayout::floorHost(const Index *Floor){
    if(not Floor->hasFloor())DEVABORT("asking for floor with non-floor index: "+Floor->strNode());
    auto iter=_floorSetupHost.find(Floor->hash());
    if(iter==_floorSetupHost.end())return Parallel::none;
    return iter->second;
}

void ParallelLayout::setFloorHosts(const Index* Idx){
    if(MPIwrapper::Size()==1)return;
    std::vector<const Index*> all;
    const Index* ix=Idx->root();
    while((ix=ix->nodeNext()))
        if(ix->hasFloor())all.push_back(ix);
    if(all.size()==0)DEVABORT("index w/o floors");

    if(_permute.size()>1){
        thisPerm=_permute;
        std::stable_sort(all.begin(),all.end(),smallerIndex);
    }

    // split into equal size()'s
    int part=int(Idx->size()/MPIwrapper::Size()+1);
    int host=0,cur=0;
    for(const Index* ix: all){
        cur+=ix->size();
        if(cur>part*(host+1)){
            host++;
        }
        _floorSetupHost[ix->hash()]=host;
    }
}

void ParallelLayout::sort(std::vector<ParallelCross *> & Cross){
    if(Cross.size()==0)return;
    //    if(Cross[0]->index()->root()->isHybrid()){
    //        PrintOutput::warning("hybrid "+Cross[0]->index()->root()->coordinates()+", ParallelLayout will not be applied");
    //        return;
    //    }
    if(Cross[0]->index()->root()->hierarchy()!=_hier)
        ABORT("Layout not for this vector: "+Cross[0]->index()->root()->hierarchy()+" != "+_hier);

    if(_permute.size()>1){
        thisPerm=_permute;
        // determine sorting on master, then broadcast to all
        std::vector<int> srt;
        std::vector<ParallelCross*>save(Cross);
        if(MPIwrapper::isMaster()){
            //CAUTION: std::sort is tricky, can fail here, found stable_sort not to show that issue
            // may be related to subtle points on the comparison function required by std::sort
            std::stable_sort(Cross.begin(),Cross.end(),ParallelLayout::smallerCross);
            for(auto c: Cross)srt.push_back(std::find(save.begin(),save.end(),c)-save.begin());
        }
        else
            srt.resize(save.size());
        MPIwrapper::Bcast(srt.data(),srt.size(),MPIwrapper::master());
        for(size_t k=0;k<srt.size();k++)Cross[k]=save[srt[k]];
    }
}
