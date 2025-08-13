// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "coefficientsLocal.h"

#include "coefficientsGlobal.h"
#include "parallel.h"
#include "parallelContinuity.h"

#include "mpiWrapper.h"
#include "timer.h"

using namespace std;

bool CoefficientsLocal::local(const Coefficients*C){
    return Parallel::owner(C->idx())==MPIwrapper::Rank() or Parallel::owner(C->idx())==Parallel::none;
}

std::complex<double> CoefficientsLocal::_dummyStorage;
static void purgeNoFloor(Coefficients* C){
    if(C->idx()->hasFloor())return;

    for(size_t k=C->childSize();k>0;k--){
        purgeNoFloor(C->child(k-1));
        if(C->child(k-1)->isLeaf() and not C->child(k-1)->idx()->hasFloor())C->childErase(k-1);
    }
}

// likely this can be fully replaced with purgeNoFloor()
static void purgeLocal(Coefficients* c, const Coefficients* orig) {
    // remove branches that do not lead to floors
    if(c->idx()->isHybrid()) { // special hybrid treatment: purge neutral and ionic subtrees individually
        // (tree algorithms won't work with non-uniform trees)
        unsigned neutralHeight = orig->child(0)->height();
        unsigned ionicHeight = orig->child(1)->height();
        if(c->childSize() > 1) { // some MPI-nodes may not have the Neutral branch
            //            c->child(0)->purge(neutralHeight);
            //            c->child(1)->purge(ionicHeight);
            purgeNoFloor(c->child(0));
            purgeNoFloor(c->child(1));
        }
        else {
            //            c->child(0)->purge(ionicHeight);
            purgeNoFloor(c->child(0));
        }
    }
    else {
        purgeNoFloor(c);
        //        c->purge(orig->height());
    }
}

CoefficientsLocal::CoefficientsLocal(const Index* I,complex<double>Val)
    :CoefficientsLocal()
{
    Coefficients C(I,Val);
    C.subTree(local,this);
    ::purgeLocal(this, &C);
    complex<double>*NextData(0);
    treeOrderStorage(NextData);
    unsetOrderedData(); // only floors are sorted according to idx()
    newCont.reset(new ParallelContinuity(this));
    _size=setSize(this);
}

int CoefficientsLocal::setSize(const Coefficients *C){
    if(C->isLeaf()){
        if(C->idx()->hasFloor())return C->idx()->size();
        // coefficient has been truncated completly
        return 0;
    }
    int _size=0;
    for(int k=0;k<C->childSize();k++)_size+=setSize(C->child(k));
    return _size;
}

CoefficientsLocal::CoefficientsLocal(const CoefficientsLocal & Other)
{
    Coefficients C(Other.idx());
    C.subTree(local,this);
    _size=Other.size();

    if(not Other.idx())DEVABORT("idx=0");
    ::purgeLocal(this, &C);
    complex<double>* NextData(0);
    treeOrderStorage(NextData);
    unsetOrderedData(); // only floors are sorted according to idx()
    newCont.reset(new ParallelContinuity(this));
    _size=setSize(this);
    operator=(Other);
    if(not idx())DEVABORT("no index");
}

CoefficientsLocal & CoefficientsLocal::operator=(const CoefficientsLocal & Other){
    if(_size!=Other._size)DEVABORT("sizes do not match");
    Coefficients::operator=(Other);
    return *this;
}

TIMER(makeContinuous,)
void CoefficientsLocal::makeContinuous(double Scal){
    STARTDEBUG(makeContinuous);
    newCont->apply(this,Scal);
    STOPDEBUG(makeContinuous);
}

/// in contiguous storage (hopefully for the same parallel layout)
CoefficientsLocal* CoefficientsLocal::cwiseProduct(const CoefficientsLocal &Rhs){
    for(int k=0;k<size();k++)anyData()[k]*=Rhs.anyData()[k];
    return this;
}

map<std::string,CoefficientsLocal*> CoefficientsLocal::_views;

CoefficientsLocal* CoefficientsLocal::view(Coefficients* C)
{
    if(_views.count(C->hash())==1){
        return _views[C->hash()];
    }
    CoefficientsLocal* v = new CoefficientsLocal();
    C->subTree(local,v,true); // view on subtree of local nodes
    ::purgeLocal(v, C);
    v->unsetOrderedData(); // only floors are Index-ordered
    v->_size=setSize(v);
    v->newCont.reset(new ParallelContinuity(v));
    _views[C->hash()]=v;
    return v;
}

std::complex<double> * CoefficientsLocal::storageData(){
    return firstLeaf()->floorData();
}

double CoefficientsLocal::norm() const{
    double nrm=const_cast<CoefficientsLocal*>(this)->Coefficients::norm();
    MPIwrapper::AllreduceMAX(&nrm,1);
    return nrm;
}

std::string CoefficientsLocal::strNode(int Precision) const{
    return Coefficients::strNode(Precision);
}
