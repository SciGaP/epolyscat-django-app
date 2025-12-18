// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "discretizationCoupledLM.h"

#include "abort.h"
#include "index.h"
#include "useMatrix.h"
#include "operatorZtimesId.h"
#include "printOutput.h"
#include "coefficients.h"
#include "purgeIndex.h"

#ifdef _USE_GSL_
#include <gsl/gsl_sf_coupling.h>
#endif


double clebgor(int l1, int l2, int L, int m1, int m2, int M){
#ifndef _USE_GSL_
    ABORT("for running TransformToCoupledLM need to compile with -D_USE_GSL_");
#else
    return std::pow(-1,l1-l2+M)*std::sqrt(2.*L+1.)*gsl_sf_coupling_3j(2*l1,2*l2,2*L,2*m1,2*m2,-2*M);
#endif
}

static std::vector<int> generateMs(int Mmax){
    std::vector<int> result;
    result.push_back(0);
    for(int i=1; i<=Mmax; i++){
        result.push_back(i);
        result.push_back(-i);
    }
    return result;
}

static std::vector<int> generateRange(int Min, int Max){
    std::vector<int> result;
    for(int i=Min; i<=Max; i++){
        result.push_back(i);
    }
    return result;
}

Index* DiscretizationCoupledLM::setupIndexStructure(int M1, int L1, int M2, int L2, std::vector<int> idx, const Index* model){
    // M level
    if(idx.size()==0){
        Index* result = new Index();
        result->setBasis(new IndexedBasisAbstract("M", generateMs(M1+M2)));
        result->setAxisName("CoupledM");

        for(int i=0; i<=(M1+M2); i++){

            std::vector<int> tmp(idx);
            tmp.push_back(i);
            result->childAdd(
                setupIndexStructure(M1, L1, M2, L2, tmp, model)
            );

            if(i!=0){
                std::vector<int> tmp2(idx);
                tmp2.push_back(-i);
                result->childAdd(
                    setupIndexStructure(M1, L1, M2, L2, tmp2, model)
                );
            }

        }
        return result;
    }else
    // L level
    if(idx.size()==1){
        Index* result = new Index();
        result->setBasis(new IndexedBasisAbstract("L", generateRange(std::abs(idx[0]), L1+L2)));
        result->setAxisName("CoupledL");

        for(int i=std::abs(idx[0]); i<=L1+L2; i++){
            std::vector<int> tmp(idx);
            tmp.push_back(i);
            result->childAdd(
                setupIndexStructure(M1, L1, M2, L2, tmp, model)
            );
        }
        return result;

    }else
    // l1 level
    if(idx.size()==2){
        Index* result = new Index();
        result->setBasis(new IndexedBasisAbstract("assocLegendre",
                generateRange(std::max(0, idx[1]-L2), std::min(L1, idx[1]+L2))));
        result->setAxisName("Eta1");

        for(int i=std::max(0,idx[1]-L2); i<=std::min(L1,idx[1]+L2); i++){
            std::vector<int> tmp(idx);
            tmp.push_back(i);
            result->childAdd(
                setupIndexStructure(M1, L1, M2, L2, tmp, model)
            );
        }
        return result;
    }else
    // l2 level
    if(idx.size()==3){
        Index* result = new Index();
        result->setBasis(new IndexedBasisAbstract("assocLegendre",
                generateRange(std::abs(idx[1]-idx[2]), std::min(L2, idx[1]+idx[2]))));
        result->setAxisName("Eta2");

        for(int i=std::abs(idx[1]-idx[2]); i<=std::min(L2,idx[1]+idx[2]); i++){
            std::vector<int> tmp(idx);
            tmp.push_back(i);
            result->childAdd(
                setupIndexStructure(M1, L1, M2, L2, tmp, model)
            );
        }
        return result;
    }
    // Below l2
    else{
        return new Index(*model);
    }
}

static bool findAngularIndices(const Index* Idx, int& m1, int& m2, int& l1, int& l2){
    m1=INT_MIN;
    m2=INT_MIN; 
    l1=INT_MIN; 
    l2=INT_MIN; 

    for(const Index* idx=Idx; idx->parent()!=0; idx=idx->parent()){
        if     (idx->parent()->axisName()=="Phi1") m1=idx->physical();
        else if(idx->parent()->axisName()=="Eta1") l1=idx->physical();
        else if(idx->parent()->axisName()=="Phi2") m2=idx->physical();
        else if(idx->parent()->axisName()=="Eta2") l2=idx->physical();
    }

    return((m1 != INT_MIN) and (m2 != INT_MIN) and (l1 != INT_MIN) and (l2 != INT_MIN));
}

static bool findCoupledIndices(const Index* Idx, int& M, int& L, int& l1, int& l2){
    M=INT_MIN;
    L=INT_MIN; 
    l1=INT_MIN; 
    l2=INT_MIN; 

    for(const Index* idx=Idx; idx->parent()!=0; idx=idx->parent()){
        if     (idx->parent()->axisName()=="CoupledM") M=idx->physical();
        else if(idx->parent()->axisName()=="Eta1") l1=idx->physical();
        else if(idx->parent()->axisName()=="CoupledL") L=idx->physical();
        else if(idx->parent()->axisName()=="Eta2") l2=idx->physical();
    }

    return((M != INT_MIN) and (L != INT_MIN) and (l1 != INT_MIN) and (l2 != INT_MIN));
}

static int findAngularDepth(const Index* Idx){
    unsigned int depth =    Idx->axisIndex("Phi1")->depth();
    depth = std::max(depth, Idx->axisIndex("Eta1")->depth());
    depth = std::max(depth, Idx->axisIndex("Phi2")->depth());
    depth = std::max(depth, Idx->axisIndex("Eta2")->depth());
    depth = depth + 1;

    if(depth != 4) PrintOutput::warning("Unsupported axes setup!");
    return depth;
}

DiscretizationCoupledLM::MapFromParent::MapFromParent(const Index* IIndex, const Index* JIndex, double* floor_factor):
    OperatorTree("Transformation to coupled angular momenta", IIndex, JIndex){

    double clebsch_gordan;

    int m1, m2, l1, l2;
    bool success = findAngularIndices(jIndex, m1, m2, l1, l2);
    int M, L, l1_, l2_;
    success = success && findCoupledIndices(iIndex, M, L, l1_, l2_);

    if(success){
        if(l1_!=l1) return;
        if(l2_!=l2) return;

        clebsch_gordan = clebgor(l1,l2,L, m1,m2,M);
        floor_factor = &clebsch_gordan;
    }


    if(floor_factor!=0){
        if(*floor_factor==0) return;

        if(iIndex->hasFloor()){
            oFloor = new OperatorZtimesId(*floor_factor,iIndex->sizeStored(), "TODO");
        }else{
            if(iIndex->childSize()!=jIndex->childSize()) ABORT("Index mismatch");
            for(int i=0; i<iIndex->childSize(); i++){
                childAdd(new MapFromParent(iIndex->child(i), jIndex->child(i), floor_factor));
            }
        }
    }else{

        for(int i=0; i<iIndex->childSize(); i++){
            for(int j=0; j<jIndex->childSize(); j++){
                OperatorTree* c = new MapFromParent(iIndex->child(i), jIndex->child(j), floor_factor);
                if(c->isZero()){
                    delete c;
                }else{
                    childAdd(c);
                }
            }
        }

    }
}

DiscretizationCoupledLM::MapToParent::MapToParent(const Index* IIndex, const Index* JIndex, double* floor_factor):
    OperatorTree("Transformation from coupled angular momenta", IIndex, JIndex){

    double clebsch_gordan;

    int m1, m2, l1, l2;
    bool success = findAngularIndices(iIndex, m1, m2, l1, l2);
    int M, L, l1_, l2_;
    success = success && findCoupledIndices(jIndex, M, L, l1_, l2_);

    if(success){
        if(l1_!=l1) return;
        if(l2_!=l2) return;

        clebsch_gordan = clebgor(l1,l2,L, m1,m2,M);
        floor_factor = &clebsch_gordan;
    }


    if(floor_factor!=0){
        if(*floor_factor==0) return;

        if(iIndex->hasFloor()){
            oFloor = new OperatorZtimesId(*floor_factor,iIndex->sizeStored(), "TODO");
        }else{
            if(iIndex->childSize()!=jIndex->childSize()) ABORT("Index mismatch");
            for(int i=0; i<iIndex->childSize(); i++){
                childAdd(new MapToParent(iIndex->child(i), jIndex->child(i), floor_factor));
            }
        }
    }else{

        for(int i=0; i<iIndex->childSize(); i++){
            for(int j=0; j<jIndex->childSize(); j++){
                OperatorTree* c = new MapToParent(iIndex->child(i), jIndex->child(j), floor_factor);
                if(c->isZero()){
                    delete c;
                }else{
                    childAdd(c);
                }
            }
        }

    }
}

DiscretizationCoupledLM::DiscretizationCoupledLM(const Discretization* Parent){
    const Index* _idx = Parent->idx();

    int depth = findAngularDepth(_idx);

    int M1 = 0;
    int M2 = 0;
    int L1 = 0;
    int L2 = 0;
    for(const Index* i=_idx->descend(depth); i!=0; i=i->nodeNext()){
        int m1, m2, l1, l2;
        findAngularIndices(i, m1, m2, l1, l2);
        M1 = std::max(M1, std::abs(m1));
        M2 = std::max(M2, std::abs(m2));
        L1 = std::max(L1, l1);
        L2 = std::max(L2, l2);
    }

    Index* innerIndex = setupIndexStructure(M1, L1, M2, L2, std::vector<int>(),_idx->descend(depth));
    innerIndex->resetFloor(_idx->firstFloor()->depth());
    innerIndex->sizeCompute();

    std::cerr<<"Setting up mapFromParent... ";
    _mapFromParent.reset(new MapFromParent(innerIndex, _idx, 0));
    std::cerr<<"Setting up mapToParent... ";
    _mapToParent.reset(new MapToParent(_idx, innerIndex, 0));

    std::cerr<<"Sizes: "<<innerIndex->sizeStored()<<"x"<<_idx->sizeStored()<<" ";
    PurgeIndex(innerIndex)
        .usingOperatorTree(dynamic_cast<const OperatorTree*>(mapFromParent()))
        .usingOperatorTree(dynamic_cast<const OperatorTree*>(mapToParent()))
        .run();
    std::cerr<<"Purged sizes: "<<innerIndex->sizeStored()<<"x"<<_idx->sizeStored()<<" ";


    std::cerr<<"Checking... ";
    Coefficients c(mapFromParent()->jIndex);
    Coefficients c1(mapFromParent()->iIndex);
    Coefficients c2(mapFromParent()->jIndex);

    c.treeOrderStorage();
    c1.treeOrderStorage();
    c2.treeOrderStorage();

    for(int i=0; i<100; i++){
        c.setToRandom();
        mapFromParent()->apply(1., c,  0., c1);
        mapToParent()->apply(1., c1, 0., c2);

        c-=c2;
        if(c.maxCoeff()>1.e-12) 
            PrintOutput::warning("WARNING! Transformation to LM not unitary: error="
                    +std::to_string(c.maxCoeff())+"\n");
    }

    std::cerr<<"done"<<std::endl;

    /*
     * Set attributes in Discretization. This is only a hint towards a full implementation
     * TODO! If this is implemented properly the index structure can be set up directly.
     */
    name = "CoupledLM("+Parent->name+")";
    _idx=innerIndex;

//    hierarchy = std::vector<std::string>();
    axis = std::vector<Axis>();
}

