// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorDiagonal.h"

#include <string>
#include "index.h"
#include "coefficients.h"
#include "coefficientsFloor.h"


using namespace std;

OperatorDiagonal::~OperatorDiagonal(){
    // for(unsigned int k=0;k<childSize();k++)delete child(k);
}

unsigned int OperatorDiagonal::vals() const{
    unsigned int n=dVal.size();
    for(size_t k=0;k<childSize();k++)n+=const_cast<OperatorDiagonal*>(this)->child(k)->vals();
    return n;
}

void OperatorDiagonal::setVal(unsigned int K, std::complex<double> Val){
    if(dVal.size()>K){
        dVal[K] = Val;
        return;
    }
    unsigned int n=0;
    for(size_t k=0;k<childSize();k++){
        int m=const_cast<OperatorDiagonal*>(this)->child(k)->vals();
        if((n+m)>K){
            child(k)->setVal(K-n, Val);
            return;
        }else
            n+=m;
    }
    ABORT(Sstr+"value of K="+K+" exceeds diagonal size at node"+index()+":"+idx()->strNode());
}

complex<double> OperatorDiagonal::val(unsigned int K) const{
    if(dVal.size()>K)return dVal[K];
    unsigned int n=0;
    for(size_t k=0;k<childSize();k++){
        int m=const_cast<OperatorDiagonal*>(this)->child(k)->vals();
        if((n+m)>K)
            return const_cast<OperatorDiagonal*>(this)->child(k)->val(K-n);
        else
            n+=m;
    }
    ABORT(Sstr+"value of K="+K+" exceeds diagonal size");
    return 0.;
}
void OperatorDiagonal::setFunctionValue(unsigned int K, std::complex<double> ValK)
{
    updateFunction(0.,doNotUpdate);
    if(dVal.size()>K){
        if(dVal.size()!=funcVal.size())funcVal=dVal;
        funcVal[K]=ValK;
        return;
    }

    unsigned int n=0;
    for(size_t k=0;k<childSize();k++){
        int m=child(k)->vals();
        if((n+m)>K){
            child(k)->setFunctionValue(K-n,ValK);
            return;
        }
        else
            n+=m;
    }
    ABORT(Sstr+"value of K="+K+"exceeds diagonal size");

}

OperatorDiagonal::OperatorDiagonal(const string & Name, const Index *Idx, timeDependentFunction Func)
    :OperatorAbstract(Name,Idx,Idx),func(Func){}

void OperatorDiagonal::add(const std::vector<complex<double> > & DVal,
                           const Index *Idx,
                           complex<double> Factor)
{
    if(iIndex->hasFloor()
            or iIndex==Idx // to work with new, multi-block structure of spectral value index
            )
    {
        if(DVal.size()!=iIndex->size()){
            ABORT(Sstr+"cannot add: data size does not match Index size: "+DVal.size()+iIndex->size());
        }
        if(dVal.size()==0)dVal.assign(DVal.size(),0.);
        for(unsigned int k=0;k<dVal.size();k++)dVal[k]+=DVal[k]*Factor;
        updateFunction(_time,keepFunction);
    }

    else
    {
        // find index in path on current operator level
        unsigned int iIns=0;
        const Index* idx=Idx;
        while(iIndex->path().size()+1<idx->path().size())idx=idx->path().back();

        // locate insert index
        while(iIndex->childSize()>iIns and idx!=iIndex->child(iIns))iIns++;
        if(iIns==iIndex->childSize())DEVABORT(Sstr+"insert index not found in operator index tree"+idx->index()+"vs"+iIndex->index());

        unsigned int iOp=0,kOp=0;
        for(;kOp<childSize();kOp++){
            // locate operator index
            if(child(kOp)->iIndex!=iIndex->child(iOp))iOp++;

            // break if found or exceeded
            if(iOp>=iIns)break;
        }

        // insert if exceeded
        if(kOp==childSize() or iOp>iIns)childAdd(new OperatorDiagonal("",idx,identityFunction));

        child(kOp)->add(DVal,Idx,Factor);
    }
}

void OperatorDiagonal::updateFunction(double Time, timeDependentFunction Func){
    if(Func!=keepFunction)func=Func;
    if(Func==doNotUpdate)return;
    for(unsigned int k=0;k<childSize();k++)child(k)->updateFunction(Time,Func);
    if(func==identityFunction)funcVal=dVal;
    else
        for(unsigned int k=0;k<dVal.size();k++)funcVal[k]=func(dVal[k],Time);
    _time=Time;
}

void OperatorDiagonal::axpy(std::complex<double> Alfa, const Coefficients &X, std::complex<double> Beta, Coefficients &Y, const double Time)
{
    if(_time!=Time){
        update(Time);
        _time=Time;
    }
    apply(Alfa,X,Beta,Y);
}

void OperatorDiagonal::apply(std::complex<double> Alfa, const Coefficients &X, std::complex<double> Beta, Coefficients &Y) const
{
    if(isLeaf()){
        if(Y.orderedData()==0 or X.orderedData()==0)DEVABORT("no orderedData() at OperatorDiagonal leaf");
        for(unsigned int k=0;k<dVal.size();k++)
            Y.orderedData()[k]=Alfa*funcVal[k]*X.orderedData()[k]+Beta*Y.orderedData()[k];
    }

    else {
        for(unsigned int k=0,l=0;k<childSize();k++,l++){
            while(child(k)->iIndex!=Y.idx()->child(l))l++;
            child(k)->apply(Alfa,*X.child(l),Beta,*Y.child(l));
        }
    }
}

string OperatorDiagonal::strNode(int Precision) const{
    string s;
    s+=tools::str(iIndex->index(),3,",")+" (diag)";    // block number
    if(childSize()!=0)
        s+=",  sub: "+tools::str(childSize());   // number of sub-blocks
    else
        if(Precision==Tree_defaultKind)
            s+=" floor("+tools::str(dVal.size())+"): "+tools::str(dVal,0,"");
        else if(Precision==Tree_ptrsSizes)
            s+=Sstr+" floor("+tools::str(dVal.size())+"): ";
        else
            s+=" floor("+tools::str(dVal.size())+"): "+tools::str(dVal,Precision,"");
    return s;
}

void OperatorDiagonal::setupAccordingToIndex(){
    if(childSize() != 0 or dVal.size() != 0) ABORT("Unsupported");

    if(iIndex->hasFloor()){
        add(std::vector<std::complex<double>>(iIndex->sizeCompute(), 0.), iIndex);
    }else{
        for(size_t i=0; i<iIndex->childSize(); i++){
            childAdd(new OperatorDiagonal("", iIndex->child(i), OperatorDiagonal::identityFunction));
            childBack()->setupAccordingToIndex();
        }
    }
}

void OperatorDiagonal::storeInCoefficients(Coefficients& c){
    if(c.idx()!=iIndex)DEVABORT("Index mismatch");
    if(childSize()==0){
        if(c.orderedData()==0)DEVABORT("structures do not match, diagonal:\n"+str()+"Coefficicients\n"+c.str());
        for(int i=0; i<vals(); i++)c.orderedData()[i] = val(i);
    }else{
        if(childSize()!=c.childSize())DEVABORT("mismatch");
        for(int i=0; i<c.childSize(); i++){
            child(i)->storeInCoefficients(*c.child(i));
        }
    }
}

void OperatorDiagonal::setFromCoefficients(const Coefficients& c){
    if(c.idx()!=iIndex)DEVABORT("Index mismatch");
    if(childSize()==0){
        if(c.orderedData()==0)DEVABORT("structures do not match, diagonal:\n"+str()+"Coefficicients\n"+c.str());
        for(size_t i=0; i<c.size(); i++){
            setVal(i,c.orderedData()[i]);
        }
    }else{
        for(int i=0; i<c.childSize(); i++){
            child(i)->setFromCoefficients(*c.child(i));
        }
    }
}
