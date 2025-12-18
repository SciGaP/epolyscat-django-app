// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "../operatorExpandIndex.h"

#include "indexNew.h"
#include "basisSub.h"
#include <set>
#include "operatorFloor.h"
#include <Eigen>
#include "eigenTools.h"

#include "indexExtract.h"
#include "indexProd.h"
#include "printOutput.h"
#include "indexAssembled.h"
#include "basisSub.h"
#include "operatorIdentity.h"

#include "parallelOperator.h"

bool OperatorExpandIndex::debug;

static bool matchLevel(const Index* Idx1, const Index* Idx2){
    if(Idx1->axisName()!=Idx2->axisName())return false;
    return (*BasisSub::superBas(Idx1->basis())==*BasisSub::superBas(Idx2->basis()));
}

OperatorExpandIndex::OperatorExpandIndex(const OperatorTree *Factor, const Index* Mdx, bool MatchI)
{
    if(Factor->idx(MatchI)->treeEquivalent(Mdx)){
        // trivial case
        nodeCopy(Factor,false);
        for(size_t k=0;k<Factor->childSize();k++)childAdd(new OperatorTree(Factor->child(k)));
        if(MatchI)replaceIndex(Mdx,new Index(*Factor->jdx()));
        else      replaceIndex(new Index(*Factor->idx()),Mdx);
        return;
    }



    // check Factor structure
    for(size_t k=1;k<Factor->childSize();k++)
        if(Factor->child(k)->idx()!=Factor->child(0)->idx() and Factor->child(k)->jdx()!=Factor->child(0)->jdx())
            DEVABORT("Algorithm assumes that one level of an operator belongs to either same row or same column");

    // local shorthand
    bool M=MatchI; // matching Index
    bool O=not M;  // other

    // IdxRight=true: Idx is to the right
    name="Id x "+name;
    if(M)iIndex=Mdx;
    else jIndex=Mdx;

    Index* odx(0);
    if(Mdx->hasFloor() and Factor->idx(M)->hasFloor() and Factor->idx(O)->hasFloor()){
        // oFloor is tensor product in the seqence as in Mdx
        if(not Factor->floor()){
            DEVABORT("arrived at product floor, but Factor is has no Floor\n"
                                        +Factor->idx(M)->str()+"\nother\n"
                                        +Factor->idx(O)->strNode()+"\nFactor\n"
                                        +Factor->str());
}
        // product of identiy and Factor floor
        size_t idSize=Mdx->size()/Factor->idx(M)->size();
        Eigen::MatrixXcd i=Eigen::MatrixXcd::Identity(idSize,idSize);
        Eigen::MatrixXcd f(Factor->floor()->matrix());

        // factor on top?
        bool facTop=matchLevel(Factor->idx(M),Mdx);
        if(facTop)oFloor=OperatorFloor::factory(std::vector<const Eigen::MatrixXcd*>({&f,&i}),name);
        else      oFloor=OperatorFloor::factory(std::vector<const Eigen::MatrixXcd*>({&i,&f}),name);

        Index* ndx=new Index();
        odx=new IndexAssembled(Factor->idx(O));
        ndx->nodeCopy(Mdx->descend(facTop),false);
        if(not facTop)std::swap(odx,ndx);
        odx->setFloor();
        for(size_t k=0;k<odx->basis()->size();k++){
            odx->childAdd(ndx->deepCopy());
        }
        delete ndx;
    }
    else{
        std::vector<int> oFun; // subset of the odx functions where branches are actually attached
        if(matchLevel(Factor->idx(M),Mdx)){
            name=Factor->name;
            odx=new IndexAssembled(Factor->idx(O));
            odx->unsetFloor();

            if(not Factor->descend())DEVABORT("factor leaf: "+Factor->str());
            if(Factor->idx(M)==Factor->descend()->idx(M)){
                // next Factor level also matches Mdx - descend
                std::vector<int> fsubs=BasisSub::subset(Factor->idx(O)->basis());
                for(size_t k=0;k<Factor->childSize();k++){
                    OperatorExpandIndex* o=new OperatorExpandIndex(Factor->child(k),Mdx,MatchI);
                    if(o->idx(O)){
                        if(o->norm()==0.)DEVABORT("zero-norm Factor block\n");
                        childAdd(o);
                        if(not k or Factor->child(k-1)->idx(O)!=Factor->child(k)->idx(O)){
                            odx->childAdd(childBack()->idx(O));
                            oFun.push_back(fsubs[Factor->child(k)->idx(O)->nSibling()]);
                        }
                    }
                    else delete o;
                }
            } else {

                // attach Idx->child(l) that contains matching function
                std::vector<size_t> matchK;
                for(size_t k=0;k<Factor->childSize();k++){
                    // function number on matching index
                    int mFun=BasisSub::subset(Factor->idx(M)->basis())[Factor->child(k)->idx(M)->nSibling()];
                    std::string mAx=Factor->idx(M)->axisName();
                    for(size_t l=0;l<Mdx->childSize();l++){
                        // add for matching function in Mdx children
                        if(mFun==BasisSub::subset(Mdx->basis())[l] and mAx==Mdx->axisName()){
                            childAdd(new OperatorExpandIndex(Factor->child(k),Mdx->child(l),MatchI));

                            // seek previous appearance of Other index
                            size_t m=0;
                            while(m<k and Factor->child(m)->idx(O)!=Factor->child(k)->idx(O))m++;
                            if(m==k){
                                // no previous - add new index child
                                if(childBack()->idx(O)){
                                    matchK.push_back(k);
                                    odx->childAdd(childBack()->idx(O));
                                    oFun.push_back(mFun);
                                }
                            }
                            else {
                                // replace with previous matching, destroy duplicate
                                const Index* d=childBack()->idx(O);
                                if(M)childBack()->replaceIndex(0,child(matchK[m])->idx(O));
                                else childBack()->replaceIndex(child(matchK[m])->idx(O),0);
                                delete d;
                            }
                            break;
                        }
                    }
                }
            }
        } else {
            // insert identity level
            name="Id x "+Factor->name;
            odx=new IndexAssembled(Mdx);
            odx->unsetFloor();

            for(size_t k=0;k<Mdx->childSize();k++){
                // add children that contain matching factor branch
                OperatorExpandIndex* o=new OperatorExpandIndex(Factor,Mdx->child(k),MatchI);
                if(o->idx(O)){
                    childAdd(o);
                    odx->childAdd(childBack()->idx(O));
                    oFun.push_back(k);
                }
                else delete o;
            }
        }
        if(not oFun.size()){
            delete odx;
            odx=0;
        } else if(oFun.size()<odx->basis()->size())
            odx->setBasis(BasisAbstract::factory(BasisSub::strDefinition(odx->basis(),oFun)));
    }

    if(not Mdx->parent() and not Factor->parent()){
        if(odx){
            odx->cleanFloors(); // remove possible duplicate floor indicators
            odx->sizeCompute();
        }
    }

    // assign the spectral (not matching) Index
    if(M)jIndex=odx;
    else iIndex=odx;

//    odx->sizeCompute();

    if(not Mdx->parent() and not Factor->parent())
        ParallelOperator::setDistribution(this);
}

Coefficients OperatorExpandIndex::coefficients(const Coefficients &Fac, const Index* IdxX){
    Coefficients C(IdxX);
    coefficients(&Fac,&C);
    if((Fac.norm()==0.)!=(C.norm()==0.))DEVABORT("bad expansion");
    return C;
}
void OperatorExpandIndex::coefficients(const Coefficients *Fac, Coefficients* C){
    if(C->isLeaf()){
        if(not Fac->isLeaf())DEVABORT("not leaf on factor");
        size_t quot=C->size()/Fac->size();
        if(matchLevel(Fac->idx(),C->idx()))
            // Fac on top - indices run slower (Eigen stores column-wise!)
            Eigen::Map<Eigen::MatrixXcd>(C->data(),quot,Fac->size())
                    =Eigen::MatrixXcd::Ones(quot,1)
                    *Eigen::Map<Eigen::MatrixXcd>(const_cast<std::complex<double>*>(Fac->data()),1,Fac->size());
        else
            // Fac on bottom - indices run faster
            Eigen::Map<Eigen::MatrixXcd>(C->data(),Fac->size(),quot)
                    =Eigen::Map<Eigen::MatrixXcd>(const_cast<std::complex<double>*>(Fac->data()),Fac->size(),1)
                    *Eigen::MatrixXcd::Ones(1,quot);
    }

    bool sameAx=matchLevel(Fac->idx(),C->idx()) and (Fac->idx()->hasFloor()==C->idx()->hasFloor());
    if(sameAx){
        std::vector<int> subF=BasisSub::subset(Fac->idx()->basis());
        for(Coefficients* cChild=C->descend();cChild!=0;cChild=cChild->nodeRight()){
            std::vector<int> subC=BasisSub::subset(cChild->parent()->idx()->basis());
            for(size_t l=0;l<Fac->childSize();l++){
                if(subC[cChild->nSibling()]==subF[l])coefficients(Fac->child(l),cChild);
            }
        }
    } else {
        for(size_t k=0;k<C->childSize();k++)
            coefficients(Fac,C->child(k));
    }
}


