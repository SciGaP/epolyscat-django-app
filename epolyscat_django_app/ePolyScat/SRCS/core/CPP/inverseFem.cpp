// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "inverseFem.h"
#include "mpiWrapper.h"
#include "index.h"
#include "coefficientsViewDeep.h"
#include "coefficientsFloor.h"
#include "operatorFloor.h"
#include "basisMat1D.h"
#include "parallelContinuity.h"
#include "basisAbstract.h"
#include "basisIntegrable.h"

using namespace std;

std::map<std::string,ParallelContinuity*> InverseFEM::_contTab;

//static bool isFem(const Index* Idx){
//    if(Idx->childSize()==1)return false;
//    const Index* idx=Idx->descend();
//    while(idx!=0 and idx->axisName()!=Idx->axisName())idx=idx->descend();
//    if(idx==0)return false;
//    return idx->basis()->integrable()!=0;
//}

InverseFEM::InverseFEM(Index* Idx, int Begin, int End)
    :Inverse("inverseFE",Idx,Idx),_correctionMap(0),_nSplit(-1),_margin(0)
{
    if(Begin<0){Begin=0;End=Idx->childSize();}

    if(Idx->hasFloor()){
        ;// do nothing
    }
    //    else if(not isFem(Idx)){
    else if(not Idx->isFem()){
        for(size_t k=0;k<Idx->childSize();k++)
            childAdd(new InverseFEM(Idx->child(k)));
    }
    else if(End-Begin==1){
        // single branch of FEM level
        childAdd(new InverseFEM(Idx->child(Begin)));
    }
    else {
        _nSplit=(Begin+End)/2;
        childAdd(new InverseFEM(Idx,Begin,_nSplit));
        childAdd(new InverseFEM(Idx,_nSplit,End));
        constructCorrection(Idx);
    }

    if(Begin==0 and size_t(End)==Idx->childSize())
        Idx->setInverseOverlap(this);
    else
        Idx->setInverseOverlap(0);
}

void InverseFEM::constructCorrection(Index *Idx){

    _margin=new Coefficients(Idx);
    _contMargin=new ParallelContinuity(_margin,_nSplit);
    _contMargin->setMargin(1.);

    Coefficients s0InvM(Idx);
    ParallelContinuity contS0(&s0InvM,_nSplit);

    // get  (T0+T1)^-1 marg
    child(0)->apply(-1.,*_margin,0.,s0InvM);
    child(1)->apply( 1.,*_margin,1.,s0InvM);
    contS0.halfDiffMargin(_contMargin); // get half the difference of margin values into _margin
    _contMargin->halfDiffMargin(_contMargin);

    // construct the map
    _correctionMap=new CorrectionMap(&s0InvM,_margin,1.);
}
InverseFEM::CorrectionMap::CorrectionMap(Coefficients *S0InvM, Coefficients *Marg, complex<double> Multi)
    :OperatorTree("correctionMap",S0InvM->idx(),Marg->idx())
{
    Eigen::MatrixXcd mOvr;
    if(Marg->root()==Marg)mOvr=Eigen::MatrixXcd::Constant(S0InvM->childSize(),Marg->childSize(),1.); // the split level - all may be connected
    //    else if(isFem(iIndex))mOvr=Eigen::MatrixXcd::Identity(S0InvM->childSize(),Marg->childSize());
    else if(iIndex->isFem())mOvr=Eigen::MatrixXcd::Identity(S0InvM->childSize(),Marg->childSize());
    else                    mOvr=BasisMat1D("1",iIndex->basis(),jIndex->basis()).mat();

    if(mOvr.size()>0){
        for(size_t i=0;i<S0InvM->childSize();i++){
            if(S0InvM->child(i)->isZero())continue;
            for(size_t j=0;j<Marg->childSize();j++){
                if(mOvr(i,j)==0.)continue;
                if(Marg->child(j)->isZero())continue;
                childAdd(new CorrectionMap(S0InvM->child(i),Marg->child(j),Multi*mOvr(i,j)));
            }
        }
    }
    if(Marg->isLeaf() or S0InvM->isLeaf()){
        S0InvM->idx()->bottomExpand();
        Marg->idx()->bottomExpand();
        if(Marg->isLeaf()!=S0InvM->isLeaf())DEVABORT("something wrong");
        floor()=new Floor(S0InvM,Marg,Multi);
//        S0InvM->idx()->bottomUnexpand();
//        Marg->idx()->bottomUnexpand();
    }

}

InverseFEM::Floor::Floor(Coefficients *S0InvM, Coefficients *RtS0invR, std::complex<double> Multi)
    :OperatorFloor(S0InvM->size(),RtS0invR->size(),"InverseFem"),
      iIndex(S0InvM->idx()),jIndex(RtS0invR->idx())
{

    // _marginDepth: where coefficients of either upper or lower margin are all non-zero
    const Index* jdx;
    for(jdx=jIndex;jdx!=0;jdx=jdx->descend())
    {
        const BasisIntegrable * bas=jdx->basis()->integrable();
        if(bas!=0){
            // on margin level, lower or upper margin data are all non-zero
            complex<double> * datBeg=RtS0invR->data()+jdx->child(bas->lowerMargin())->posIndex(jIndex);
            complex<double> * datEnd=datBeg          +jdx->child(bas->lowerMargin())->sizeStored();
            if(std::find(datBeg,datEnd,0.)==datEnd){_upperMargin=false;break;}

            datBeg=RtS0invR->data()+jdx->child(bas->upperMargin())->posIndex(jIndex);
            datEnd=datBeg          +jdx->child(bas->upperMargin())->sizeStored();
            if(std::find(datBeg,datEnd,0.)==datEnd){_upperMargin=true;break;}
        }
    }
    if(jdx==0)DEVABORT("failed to determine _marginDepth - check algorithm");
    _marginDepth=jdx->depth();

    // get zinv= [R^T S0^-1 R]^-1 (on margin only)
    vector<complex<double> > zinv(jIndex->sizeCompute(),0.);
    for(size_t k=0;k<jIndex->sizeStored();k++)
        if(RtS0invR->data()[k]!=0.)
            zinv.data()[k]=1./RtS0invR->data()[k];

    // get -(S0^-1 R Zinv)
    _s0InvM.assign(S0InvM->data(),S0InvM->data()+S0InvM->size());
    vector<complex<double> > aux(_s0InvM.size(),0.);
    axpy(-0.5*Multi,zinv.data(),zinv.size(),0.,aux.data(),aux.size());
    std::swap(aux,_s0InvM);
}

static const complex<double> *pC,*pX,*pY;
void InverseFEM::Floor::axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
                             const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const
{
    if(Beta!=1.)
        for(size_t k=0;k<SizY;k++)Y[k]*=Beta;
    const complex<double>* C=_s0InvM.data();
    pC=C;pX=X;pY=Y;
    axpyRecursive(iIndex,jIndex,Alfa,C,X,Y);
}


void InverseFEM::Floor::axpyRecursive(const Index *IIdx, const Index *JIdx,const std::complex<double> & Alfa,
                                      const std::complex<double> *&C,const std::complex<double> *&X,std::complex<double> *&Y) const {
    // Algorithm (slightly tricky):
    // tree indices indices (i[0],i[1],...,i[D-1])
    // coefficient storage is contiguous and "row-wise" wrt to their multi-indices, i.e. i0 runs slowest
    // let m=marginDepth be the level of the margin
    // for i[m]=s, margin values are contiguous for sub-indices (i[m+1]...i[D-1])
    // on m, X is set to the beginning of the stretch of contiguous margin values
    // C and Y run through the same range and are incremented throughout
    // Note: C may be quite sparse, such that one may skip the sparse ranges.

    bool iBot=IIdx->isBottom();
    bool jBot=JIdx->isBottom();

    if(JIdx->depth()<size_t(_marginDepth)){
        // if the following never fails, the distinction of IIdx and JIdx can be removed
        if(IIdx->childSize()!=JIdx->childSize())
            DEVABORT("algorithm error - all floors connected to same margin must be equivalent\n"+IIdx->str()+"\nother\n"+JIdx->str());
        for(size_t k=0;k<IIdx->childSize();k++){
            if(iBot or jBot)DEVABORT("cannot use recursion at bottom");
            axpyRecursive(IIdx->child(k),JIdx->child(k),Alfa,C,X,Y);
        }
    }
    else {
        // the following allows non-product floor indices
        const complex<double> * margX=X; // beginning of contiguous part of margin
        if(_upperMargin)for(size_t k=0;k<JIdx->basis()->integrable()->upperMargin();k++)margX+=jBot?1:JIdx->child(k)->size();
        else            for(size_t k=0;k<JIdx->basis()->integrable()->lowerMargin();k++)margX+=jBot?1:JIdx->child(k)->size();

        // below we should keep track of sparsity and reduce the operation correspondingly
        for(size_t k=0;k<IIdx->childSize();k++)
            for(size_t l=0;l<(iBot?1:IIdx->child(k)->size());l++,Y++,C++){
                *Y+=*C*margX[l]*Alfa;
            }
        X+=JIdx->size(); // increment to next segment of rhs vector
    }
}

void InverseFEM::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{

    // level without continuity below (will be floor in general)
    if(childSize()==0)
        iIndex->localInvOvr()->apply(A,Vec,B,Y);

    // standard tree level (but above
    else if(_nSplit==-1)
        for(size_t k=0;k<childSize();k++)child(k)->apply(A,*Vec.child(child(k)->jIndex->nSibling()),B,*Y.child(child(k)->iIndex->nSibling()));

    // split FE level
    else{
        for(size_t k=0;k<childSize();k++)child(k)->apply(A,Vec,B,Y);
        applyCorrection(1.,Vec,B,Y);
    }
}

void InverseFEM::applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
{
    if(_correctionMap==0){
        ABORT("remains to be implemented");
        _contTab[Y.hash()]->apply(&Y,1.);
    }
    else {
        // get the split view on the Y
        ParallelContinuity *_contY=_contTab[Y.hash()+"B"+tools::str(_nSplit)];
        if(_contY==0){
            // set up if needed
            _contY=new ParallelContinuity(&Y,_nSplit);
            _contTab[Y.hash()+"B"+tools::str(_nSplit)]=_contY;
        }
        // extract difference of margin values into _margin and apply correction
        _contY->halfDiffMargin(const_cast<ParallelContinuity*>(_contMargin));
        _correctionMap->apply(1.,*_margin,1.,Y);
    }
}
