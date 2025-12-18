// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
//#include "operator.h"

#include "printOutput.h"

//resolve forward declarations
#include "operatorSingle.h"
#include "coefficients.h"
#include "coefficientsFloor.h"
#include "discretization.h"
#include "index.h"
using namespace std;
using namespace tools;
//#include "eigenNames.h"


deque<UseMatrix> OperatorSingle::matsTable(0); /// table to contain all mats of operator single
map<string,UseMatrix> OperatorSingle::matsTableNew; /// table to contain all mats of operator single

const UseMatrix* OperatorSingle::matsAdd(UseMatrix &mat, string Hash){
    //check if an identical matrix has already been added
    if(Hash=="noHash"){
        for(unsigned int k=0;k<matsTable.size();k++){
            if(matsTable[k].cols()==0 or matsTable[k].rows()==0) ABORT("zero size matrix");
            if(matsTable[k]==mat)return (&(matsTable[k]));
        }
        //did not find identical matrix
        matsTable.push_back(mat);
        return (&(matsTable.back()));
    }
    else {
        if(matsTableNew.count(Hash)==0)matsTableNew[Hash]=mat;
        if(not (matsTableNew[Hash]==mat)){
            if( not (mat-matsTableNew[Hash]).isZero(matsTableNew[Hash].maxAbsVal()*1.e-12)){
                PrintOutput::warning("non-unique Hash: "+Hash,10);
                unsigned int k=0;
                for(;matsTableNew.count(Hash+tools::str(k))!=0;k++);
                Hash+=tools::str(k);
                matsTableNew[Hash]=mat;
            }
        }
        return &matsTableNew[Hash];
    }
}

void OperatorSingle::fuse(std::vector<OperatorSingle*> &Ops){
    if(Ops.size()<2)return;
    for (unsigned int k=0;k<Ops.size();k++){
        if(Ops[k]==0)ABORT("this should not happen");
        for(unsigned int l=Ops.size()-1;l>k;l--){
            if(Ops[k]->absorb(Ops[l]))Ops.erase(Ops.begin()+l);
        }
    }
}

OperatorSingle::~OperatorSingle(){
}

void OperatorSingle::apply(std::complex<double>  * InOut) const //Y = Operator*X+const;
{
    ABORT("not default implementation for apply");
}

void OperatorSingle::apply(std::vector<std::complex<double> > & InOut) const //Y = Operator*X+const;
{
    ABORT("not default implementation for apply");
}

const Discretization* OperatorSingle::dataDisc(const Discretization *idisc, const Discretization *jdisc){
    if(idisc==jdisc)return idisc;
    if(idisc->parent==jdisc)return idisc;
    if(jdisc->parent==idisc)return jdisc;
    ABORT("OperatorSingle::dataDisc: cannot assign data: two discretizations not in parent-child relation");
    return idisc;
}

double OperatorSingle::norm() const{
    if(oNorm<0.)
        const_cast<OperatorSingle*>(this)->oNorm
            =const_cast<OperatorSingle*>(this)->setNorm();
    return oNorm;
}

double OperatorSingle::setNorm() {
    oNonzeros=0;
    oNorm=1.;
    for (unsigned int k=0;k<mats.size();k++){
        oNorm*=mats[k]->maxAbsVal();
        oNonzeros+=mats[k]->nonZeros();
    }
    if(mats.size()==0)oNorm=0.;
    return oNorm;
}

OperatorSingle::OperatorSingle(const string Name, const string Def,  const Index * IIndex, const Index * JIndex):
    name(Name),definition(Def),iIndex(IIndex),jIndex(JIndex),oNorm(-1.),timeDepFac(0){}

string OperatorSingle::strStructure() const{
    string s=name+": "+definition+" ("
            +tools::str(iIndex->index())+"|"
            +tools::str(jIndex->index())+")";
    return s;
}

void OperatorSingle::matrix(UseMatrix & mat) const{
    DEVABORT("needs re-implementation");
    // default matrix setup through axpy operation
    mat=UseMatrix::Zero(iIndex->sizeStored(),jIndex->sizeStored());
    vector<complex<double> > xS(jIndex->sizeStored(),0.),yS(iIndex->sizeStored(),0.);
    CoefficientsFloor x(jIndex,xS.data()),y(iIndex,yS.data());
    // directly write - full matrix storage
    for(unsigned int j=0,j0=0;j<jIndex->sizeStored();j++,j0+=iIndex->sizeStored()){
        xS.assign(xS.size(),0.);
        yS.assign(yS.size(),0.);
        xS.data()[j]=1.;
        axpy(x,y);
        for (unsigned int i=0;i<iIndex->sizeStored();i++) {mat.data()[j0+i]=yS.data()[i];}
    }
}
