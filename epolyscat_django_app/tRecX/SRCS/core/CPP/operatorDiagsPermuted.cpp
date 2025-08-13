// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "operatordDiagsPermuted.h"

#include "operatorTree.h"
#include "operatorRALL.h"
#include "operatorPermute.h"
#include "operatorFloor.h"
#include "index.h"
#include "qtEigenDense.h"
#include "basisVector.h"
#include <Sparse>

#include "operatorRALL.h"

#include "eigenTools.h"

OperatorDiagsPermuted::~OperatorDiagsPermuted(){
    delete perm;
    delete back;
    delete diag;
}

OperatorDiagsPermuted::OperatorDiagsPermuted(OperatorTree& Op, double DiagonalThreshold)
    :OperatorTree(Op.name,Op.iIndex,Op.jIndex),perm(0),back(0),diag(0)
{
    // create list of diagonals on child-level, break if non-diagonal block
    int iDepth=Op.iIndex->depth(),jDepth=Op.jIndex->depth();
    std::vector<OperatorTree*> blocks;
    std::map<const Index*,std::map<const Index*,Eigen::VectorXcd>> diags;

    for(OperatorTree* op=&Op;
        op!=0
        and op->iIndex->depth()-iDepth<2 //DEBUG only(?): what is that number 2 here doing?
        and op->jIndex->depth()-jDepth<2;
        op=op->nodeNext()){
        if(op->iIndex->depth()-iDepth and op->jIndex->depth()-jDepth){
            Eigen::MatrixXcd m=op->matrix();
            if(not m.isDiagonal(DiagonalThreshold)){diags.clear();break;}
            diags[op->iIndex][op->jIndex]=m.diagonal();
        }
    }
    if(diags.size())construct(diags);
}

OperatorDiagsPermuted::OperatorDiagsPermuted(std::string Name, const Index *IIndex, const Index*JIndex,
                                             std::map<const Index*,std::map<const Index*,Eigen::VectorXcd>> &Diags)
    :OperatorTree(Name,IIndex,JIndex),perm(0),back(0),diag(0)
{
    construct(Diags);
}

static void placeInto(Index* Idx, Index* IAdd, const std::vector<unsigned int> & Place){
    Index* idx=Idx->nodeAt({Place.begin(),Place.end()-1});
    for(size_t k=idx->childSize();k<Place.back()+1;k++)idx->childAdd(new Index());
    idx->childReplace(Place.back(),IAdd);
}

void OperatorDiagsPermuted::construct(std::map<const Index *, std::map<const Index *, Eigen::VectorXcd> > &Diags)
{
    if(not Diags.size())DEVABORT("no diags");
    if(not Diags.begin()->second.size())DEVABORT("empty row");

    while(Diags.size()){
        const Index* iParent=Diags.begin()->first->parent();
        const Index* jParent=Diags.begin()->second.begin()->first->parent();
        if(not iParent)DEVABORT("no iParent");
        if(not jParent)DEVABORT("no jParent");
        if(iIndex==iParent and jIndex==jParent)break;


        if(not diag)diag=new OperatorTree("fullOp",iIndex->deepCopy(iParent->depth()),jIndex->deepCopy(jParent->depth()));
        if(not perm)perm=new OperatorTree("fullPerm",diag->jdx(),jdx());
        if(not back)back=new OperatorTree("fullBack",idx(),diag->idx());

        std::map<const Index *, std::map<const Index *, Eigen::VectorXcd> >sub;
        for(size_t i=0;i<iParent->childSize();i++){
            for(size_t j=0;j<jParent->childSize();j++)
                if(Diags.count(iParent->child(i)) and Diags[iParent->child(i)].count(jParent->child(j)))
                {
                    sub[iParent->child(i)][jParent->child(j)]=Diags[iParent->child(i)][jParent->child(j)];
                    Diags[iParent->child(i)].erase(jParent->child(j));
                }
            if(not Diags[iParent->child(i)].size())Diags.erase(iParent->child(i));


        }
        OperatorDiagsPermuted* opd=new OperatorDiagsPermuted("sub",iParent,jParent,sub);
        addColumnwise(opd);
        // add to indices
        placeInto(diag->idx(),opd->diag->idx(),opd->idx()->index());
        placeInto(diag->jdx(),opd->diag->jdx(),opd->jdx()->index());

        // add top view of perm,back,diag
        perm->addColumnwise(opd->perm);
        back->addColumnwise(opd->back);
        diag->addColumnwise(opd->diag);
        opd->perm=0;
        opd->back=0;
        opd->diag=0;

        if(not Diags.size()){
            diag->idx()->sizeCompute();
            diag->jdx()->sizeCompute();
            _lhs.reset(new Coefficients(diag->idx()));
            _rhs.reset(new Coefficients(diag->jdx()));
            return;
        }
    }
    // all Diags here belong to the parent pair (iIndex,jIndex)

    perm=new OperatorPermute("perm",1,jIndex,true);
    back=new OperatorPermute("back",1,iIndex,false);
    _rhs.reset(new Coefficients(perm->iIndex));
    _lhs.reset(new Coefficients(back->jIndex));
    if(perm->iIndex->childSize()!=back->jIndex->childSize())
        DEVABORT("number of blocks in lhs and rhs does not match");


    // construct the matrices in the transposed index
    std::vector<Eigen::MatrixXcd> mat;
    for(size_t k=0;k<perm->iIndex->childSize();k++)
        mat.push_back(Eigen::MatrixXcd::Zero(back->jIndex->child(k)->size(),perm->iIndex->child(k)->size()));
    for(auto &row: Diags)
        for(auto &block: row.second){
            Eigen::VectorXcd &diag=block.second;
            for(size_t k=0;k<mat.size();k++){
                mat[k](row.first->nSibling(),block.first->nSibling())=diag[k];
            }
        }
    diag=(new OperatorTree("diag",back->jIndex,perm->iIndex));
//    for(size_t k=0;k<perm->iIndex->childSize();k++)
//        diag->childAdd(new OperatorTree("diagblock","n/a",back->jIndex->child(k),perm->iIndex->child(k),
//                                        OperatorFloor::factory({&mat[k]},"DiagBlock")));
}

void OperatorDiagsPermuted::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    if(perm and diag and back){
        perm->apply(A,Vec,0.,*_rhs);
        diag->apply(1,*_rhs,0.,*_lhs);
        back->apply(1,*_lhs,B,Y);
    }
    else
        OperatorTree::apply(A,Vec,B,Y);
}

bool OperatorDiagsPermuted::symmetrize(){
    if(not idx()->treeEquivalent(jdx()))return false;
    if(not diag->isComplexSymmetric())return false;
    Coefficients jC(jdx()),iC(idx()),iT(diag->idx()),jT(diag->jdx());

    for(size_t k=0;k<5;k++){
        jC.setToRandom();
        perm->apply(1.,jC,0.,jT);
        iT=jT;
        back->apply(1.,iT,0.,iC);
        if(not (iC-=jC).isZero())return false;
    }

    const Index* diagIdx=diag->idx();
    diag->replaceIndex(perm->idx(),perm->idx());
    back->replaceIndex(0,perm->idx());
    if(diagIdx!=diag->idx())delete diagIdx;

    _lhs.reset(new Coefficients(diag->idx()));
    _rhs.reset(new Coefficients(diag->jdx()));

    return true;
}




