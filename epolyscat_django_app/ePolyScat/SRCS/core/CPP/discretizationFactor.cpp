// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "discretizationFactor.h"

#include <vector>
#include <memory>

#include "str.h"
#include "basisAbstract.h"
#include "basicDisc.h"
#include "coefficients.h"
#include "coefficientsFloor.h"
#include "indexNew.h"
#include "operatorTree.h"
#include "timer.h"
#include "axisTree.h"

// if top-axis is hybrid, there is an extra place-holder "Hybrid" axis for the Tree's root
// this is removed when the whole tree gets attached as a subtree
void addHybrid(AxisTree* AxTree, AxisTree* Sub){
    if(Sub->name=="Hybrid"){
        for(int k=0;k<Sub->childSize();k++){
            AxTree->childAdd(Sub->child(k));
            Sub->child(k)=0; // hide from deletion
        }
        delete Sub;
    }
    else
        AxTree->childAdd(Sub);
}

using namespace std;

DiscretizationFactor::DiscretizationFactor(const Discretization *Disc, std::string Axes, bool Complement, int noChannels):
    tempParentC(0),tempViewC(0),tempC1(0),resIdx(0),helper_ion(0),hasChannelLevel(noChannels>0),factorization(Residual_Ion)
{
    parent = Disc;

    if(noChannels>0){
        // add channels level
        _axisTree.reset(new AxisTree(Axis("Vec",noChannels,0.,double(noChannels-1),"automatic",noChannels)));
        if(Complement)addHybrid(_axisTree.get(),Disc->axisTree()->complement(tools::splitString(Axes,' '))); // add, possibly remove auxilicar "Hybrid" axis
        else          addHybrid(_axisTree.get(),Disc->axisTree()->factor(tools::splitString(Axes,' ')));
    }
    else {
        if(Complement)_axisTree.reset(Disc->axisTree()->complement(tools::splitString(Axes,' ')));
        else          _axisTree.reset(Disc->axisTree()->factor(tools::splitString(Axes,' ')));
    }
    // determine factor and quotient axes
    for(unsigned int k=0;k<Disc->getAxis().size();k++){

        // special case
        if(Disc->getAxis()[k].basDef[0].funcs.find("Lshape")!=string::npos){
            DEVABORT("Lshape in axes no longer supported");
        }
    }

    idx()=new IndexNew(_axisTree.get());
    dynamic_cast<IndexNew*>(idx())->buildOverlap();

    if(idx()->hierarchy()=="Vec.Ion")idx()->resetFloor(1);   // haCC only

    // save the views for non-complemented one - arbitrary choice
    // Create a view of the parent disc to have factor indices on top
    if(Complement) initialize_helpers(Disc);
}

DiscretizationFactor::~DiscretizationFactor()
{
    if(tempParentC!=0) {delete tempParentC; tempParentC=0;}
    if(helper_ion!=0)  {delete helper_ion;  helper_ion=0;}
}

void DiscretizationFactor::initialize_helpers(const Discretization *Disc)
{
    if(Disc->idx()->depth()!=0)ABORT("discretization index must be top level (not, e.g., block of hybrid)");

    // Create a View with channel indices on top and ion indices below
    // temporary storage
    tempParentC = new Coefficients(Disc->idx());

    // channel indices in original discretization
    vector<unsigned int > chanInd;
    for(Index* s=idx(); not s->isLeaf(); s=s->child(0)){
        if((s->axisName()=="Vec" or s->axisName()=="Channel") and s->depth()==0 and hasChannelLevel) continue;
        const Index* p=Disc->idx();
        for(; not p->isLeaf(); p=p->child(0))
            if(p->axisName()==s->axisName() and std::find(chanInd.begin(),chanInd.end(),p->depth())==chanInd.end()) break;
        chanInd.push_back(p->depth());
    }

    // Ionic indices in original discretization
    vector<unsigned int > IonInd;
    for(const Index* s=Disc->idx(); not s->isLeaf(); s=s->child(0))
        if(std::find(chanInd.begin(),chanInd.end(),s->depth())==chanInd.end())IonInd.push_back(s->depth());

    // Choose sorting
    vector<unsigned int> perm = chanInd;
    if(chanInd.back()==Disc->idx()->height()-1)    {factorization = Ion_Residual;  perm.insert(perm.begin(),IonInd.begin(),IonInd.end());}
    else if(IonInd.back()==Disc->idx()->height()-1){factorization = Residual_Ion;  perm.insert(perm.end(),IonInd.begin(),IonInd.end());}
    else                                            ABORT("Unknown case");

    // delete from bottom if no permute below original floor
    for(int k=Disc->idx()->height()-1;k>=Disc->idx()->heightAboveFloor();k--) {
        if(k!=perm[k])  break;
        perm.pop_back();
    }

    resIdx = new Index(*Disc->idx());
    tempC1 = new Coefficients(resIdx,tempParentC);

    // Permuted Coefficient
    tempViewC = new Coefficients();
    tempC1->permute(perm,*tempViewC,true);
}

void DiscretizationFactor::contractEachChannel(const Coefficients &C, const Coefficients &Ion, Coefficients &Cofactor)
{
    if(C.isZero() or Ion.isZero()) {
        Cofactor.setToZero();
        return;
    }

    switch(factorization){

    case Residual_Ion:

        if(not Cofactor.idx()->hasFloor()){
            if(C.idx()->axisName()!=Cofactor.idx()->axisName()) ABORT("Axes mis match");
            //        if(C.childSize()!=Cofactor.childSize()) ABORT("Child sizes do not match");

            for(int k=0;k<Cofactor.childSize();k++)
                contractEachChannel(*C.child(k),Ion,*Cofactor.child(k));
        }
        else{
            for(unsigned int k=0; k<Cofactor.idx()->sizeStored();k++){
                if(Ion.idx()->hierarchy().find("Phi")!=string::npos and Ion.idx()->hierarchy().find("Eta")!=string::npos
                        and Ion.idx()->hierarchy().find("Rn")!=string::npos)
                    Cofactor.floorData()[k] = conj(computeInnerProduct(*C.child(k),Ion));
                else
                    Cofactor.floorData()[k]=conj(C.child(k)->innerProductUnscaled(&Ion));
            }
        }
        return;

    case Ion_Residual:

        if(not Ion.idx()->hasFloor()){
            if(C.idx()->axisName()!=Ion.idx()->axisName()) ABORT("Axes mis match");
            if(C.childSize()!=Ion.childSize()) ABORT("Child sizes do not match");
            for(int k=0;k<Ion.childSize();k++)
                contractEachChannel(*C.child(k),*Ion.child(k),Cofactor);
        }
        else{
            for(unsigned int k=0; k<Ion.idx()->size();k++){
                Cofactor.axpy(conj(Ion.floorData()[k]),C.child(k));
            }
        }
        return;

    default:
        ABORT("Unknown case");
    }
}

complex<double> DiscretizationFactor::computeInnerProduct(const Coefficients &C, const Coefficients &Ion)
{

    if(C.idx()->basis()->isAbsorptive())return 0;

    if (C.hasFloorData()) {
        complex<double> result = C.floorInnerProduct(&Ion,false);
        return result;
    }
    else{
        complex< double > result = 0.0;

        if(C.childSize()==Ion.childSize()){
            for (unsigned int k=0; k<C.childSize(); k++){
                result += computeInnerProduct(*C.child(k),*Ion.child(k));
            }
        }
        else{  // Some restrictions were imposed

            // M restriction
            if(C.idx()->axisName().find("Phi")!=string::npos){
                if(C.idx()->axisName()==Ion.idx()->axisName() and C.childSize()==1 and C.childSize()!=Ion.childSize()){
                    DEVABORT("Index::basisSet not longer available");
                    int m=0;//C.idx()->basisSet()->parameters()[0];  // should be size 1
                    bool found = false;
                    for(unsigned int j=0;j<Ion.childSize();j++){
                        DEVABORT("Index::basisSet no longer available");
//                        if(m==Ion.idx()->basisSet()->parameters()[j]){
//                            found=true;
//                            result += computeInnerProduct(*C.child(0),*Ion.child(j));
//                            break;
//                        }
                    }
                    if(found==false) ABORT("Could not find m index in the Ion.");
                }
                else ABORT("Unknown M restriction");
            }

            // l restriction
            if(C.idx()->axisName().find("Eta")!=string::npos){
                for (unsigned int k=0; k<Ion.childSize(); k++){
                    result += computeInnerProduct(*C.child(k),*Ion.child(k));
                }
            }
        }

        return result;
    }
}

void DiscretizationFactor::contract(const Coefficients &C, const Coefficients &Ion, Coefficients &Cofactor)
{  
    if(tempParentC==0) ABORT("Indices do not match. Try calling this function with the cofactor!");
    if(tempParentC->idx()!=C.idx()) ABORT("Indices do not match. Check the order of arguments!");
    *tempParentC = C; // make a copy to be accessed by view

    // axpy with overlap, make continuous
    if(helper_ion==0) helper_ion = new Coefficients(Ion.idx());
    helper_ion->setToZero();
    //    if(Ion.idx()->hierarchy()!="Vec.Ion.NONE"){
    if(Ion.idx()->hierarchy()!="Vec.Ion"){
        Ion.idx()->localOverlap()->apply(1,Ion,1.,*helper_ion);
        helper_ion->makeContinuous();
    }
    else
        *helper_ion=Ion;   // Assumes Identity overlap. haCC only

    Cofactor.setToZero();

    if(not hasChannelLevel)
        contractEachChannel(*tempViewC,*helper_ion,Cofactor);
    else{
        if(Ion.childSize()!=Cofactor.childSize()) ABORT("Channel sizes do not match");
        for(unsigned int k=0;k<Ion.childSize();k++) contractEachChannel(*tempViewC,*helper_ion->child(k),*Cofactor.child(k));
    }

}
