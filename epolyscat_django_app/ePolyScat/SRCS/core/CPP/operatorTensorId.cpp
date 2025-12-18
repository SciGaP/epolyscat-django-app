// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorTensorId.h"

#include "basisAbstract.h"
#include "indexQuot.h"
#include "indexProd.h"
#include "operatorZGxZD.h"
#include "operatorZDxZG.h"
#include "operatorTree.h"
#include "coefficients.h"
#include "coefficientsPermute.h"

using namespace std;
OperatorTensorId::OperatorTensorId(const OperatorAbstract *Factor, const Index *JFull, const Index* IFull)
    :OperatorAbstract("Id (x) "+Factor->name,IFull,JFull),factor(Factor)
{
    // get the quotient index
    IndexQuot * iQuot = new IndexQuot(JFull,Factor->jIndex);


    // get iQuot (x) iFactor and iQuot (x) jFactor
    IndexProd* iProd=new IndexProd(iQuot,Factor->iIndex);
    IndexProd* jProd=new IndexProd(iQuot,Factor->jIndex);

    // get the permutations relative to the JFull, and, if present, IFull
    vector<unsigned int> iPerm,jPerm;
    for(int k=0;k<JFull->height();k++)
        jPerm.push_back(jProd->descend(k)->axisLevel(JFull));

    for(int k=0;IFull!=0 and k<IFull->height();k++)
        iPerm.push_back(iProd->descend(k)->axisLevel(IFull));

    iIndex=iProd;
    if(IFull!=0)iIndex=IFull;
    jIndex=JFull;
    iIndex->sizeCompute();
    jIndex->sizeCompute();

    iLoc=new CoefficientsPermute(iIndex,iPerm);
    jLoc=new CoefficientsPermute(jIndex,jPerm);

    // vectors of nodes on factor level
    Coefficients* iNode=iLoc->descend(iQuot->height());
    Coefficients* jNode=jLoc->descend(iQuot->height());
    while(iNode!=0){
        indexReplace(iNode,factor->iIndex);
        indexReplace(jNode,factor->jIndex);
        iVfac.push_back(iNode);
        jVfac.push_back(jNode);
        iNode=iNode->nodeRight();
        jNode=jNode->nodeRight();
    }

    delete iQuot;
}

void OperatorTensorId::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    // copy vector into view of factor-vector
    const_cast<CoefficientsPermute*>(jLoc)->fromOrig(Vec); // copy to factor storage

    // apply below the factor level
    for(int k=0;k<jVfac.size();k++)
        factor->apply(A,*jVfac[k],B,*iVfac[k]);

    // put into output
    iLoc->toOrig(Y); // get back from factor storage
}

void OperatorTensorId::indexReplace(Coefficients *C, const Index *I){
    if(not (*C->idx()->basis()==*I->basis())){
        cout<<"permuted\n"<<C->idx()->str()<<endl;
        cout<<"operator\n"<<I->str()<<endl;
        ABORT("indices not equivalent - cannot replace");
    }
        C->setIdx(I);
        for(int k=0;k<C->childSize();k++)
            indexReplace(C->child(k),I->child(k));
}
