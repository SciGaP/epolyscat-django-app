// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "inverseDVRmixFE.h"

#include "mpiWrapper.h"
#include "parallelContinuity.h"

#include "operatorTree.h"
#include "index.h"
#include "operatorFloor.h"
#include "printOutput.h"

using namespace std;

static void neighbors(const Index* Idx, const Index* &LoNeig, const Index* &UpNeig ){
    // this should go into a Index or Tree function
    LoNeig=0;
    UpNeig=0;
    // ascend until matching coordinate (return 0 if none), keep track of path
    string ax=Idx->axisName();
    const Index* fem=Idx;
    vector<unsigned int> path;
    while(fem->parent()!=0 and ax!=fem->parent()->axisName()){
        path.insert(path.begin(),fem->nSibling());
        fem=fem->parent();
    }
    // descend at neighbor
    if(fem!=0){
        if(fem->nSibling()>0)     LoNeig=fem->parent()->child(fem->nSibling()-1)->nodeAt(path);
        if(fem->rightSibling()!=0)UpNeig=fem->rightSibling()->nodeAt(path);
    }
}


InverseDVRmixFE::InverseDVRmixFE(Index* Idx):Inverse("DVRmixFE",Idx,Idx)
{
    Idx->setInverseOverlap(this);

    if(not Idx->hasFloor()){
        // recursively build raw inverse tree
        inv0=new OperatorTree("inv",Idx,Idx);
        for (size_t k=0;k<Idx->childSize();k++){
            InverseDVRmixFE* t=new InverseDVRmixFE(Idx->child(k));
            childAdd(t);
            inv0->childAdd(t->inv0);
        }
    }

    else
    {
        // loop through directions
        vector<Eigen::MatrixXcd> iMats;
        const Index *ldx,*udx;
        for(size_t d=0;d<=Idx->heightAboveBottom();d++){

            // find lower and upper neighbors at d'th level
            neighbors(Idx->descend(d),ldx,udx);
            for(int up=0;up<int(d);up++){if(ldx!=0)ldx=ldx->parent();if(udx!=0)udx=udx->parent();}

            Eigen::MatrixXcd loNeig,upNeig;
            if(ldx!=0)loNeig=ldx->localOverlap()->floor()->matrixFactor(d);
            if(udx!=0)upNeig=udx->localOverlap()->floor()->matrixFactor(d);

            iMats.push_back(Idx->localOverlap()->floor()->matrixFactor(d));

            //the Eigen::...Corner... function appear not to be working
            int last=iMats[d].cols()-1,loLast=loNeig.cols()-1;

            // present overlap not diagonal
            if(not iMats[d].isDiagonal()){
                if((ldx!=0 and not loNeig.isDiagonal()) or
                    (udx!=0 and not upNeig.isDiagonal())){
                    PrintOutput::warning("cannot have two subsequent FE intervals on axis "+Idx->axisName(),10,0,
                                         "...at present, only isolated FE intervals can be used that are embedded in a DVR axis");
                        ABORT("only DVR with isolated FE intervals sallowed");
                }
                // scale margins columns and rows, and add neighbors' corners
                if(ldx!=0){
                    iMats[d].leftCols(1)*=sqrt(0.5);
                    iMats[d].topRows(1)*=sqrt(0.5);
                    iMats[d](0,0)+=loNeig(loLast,loLast)*0.5;
                }
                if(udx!=0){
                    iMats[d].rightCols(1)*=sqrt(0.5);
                    iMats[d].bottomRows(1)*=sqrt(0.5);
                    iMats[d](last,last)+=upNeig(0,0)*0.5;
                }

                // invert
                iMats[d]=iMats[d].inverse();

                // multiply margins (as these will be distributed by continuity)
                if(ldx!=0){
                    iMats[d].leftCols(1)*=sqrt(2.);
                    iMats[d].topRows(1)*=sqrt(2.);
                }
                if(udx!=0){
                    iMats[d].rightCols(1)*=sqrt(2.);
                    iMats[d].bottomRows(1)*=sqrt(2.);
                }
            }

            // present overlap diagonal:
            else {
                // at diagonal neighbor's side, average
                if(ldx!=0 and loNeig.isDiagonal()){iMats[d](0,0)=0.5*(iMats[d](0,0)+loNeig(loLast,loLast));}
                if(udx!=0 and upNeig.isDiagonal()){iMats[d](last,last)=0.5*(iMats[d](last,last)+upNeig(0,0));}

                // invert
                iMats[d]=iMats[d].inverse();

                // at full neighbor's side, set to zero
                if(ldx!=0 and not loNeig.isDiagonal())iMats[d](0,0)=0.;
                if(udx!=0 and not upNeig.isDiagonal())iMats[d](last,last)=0.;
            }
        }
        // both directions done - construct OperatorFloor
        inv0=new OperatorTree("inverse",iIndex,jIndex);
        vector<const Eigen::MatrixXcd*> pMats;
        for(size_t k=0;k<iMats.size();k++)pMats.push_back(&iMats[k]);
        inv0->floor()=OperatorFloor::factory(pMats,"invers"+iIndex->hash()+jIndex->hash());
    }
}

void InverseDVRmixFE::parallelSetup() const {
    if(MPIwrapper::Size()>1)ABORT("cannot be run in parallel");
}

void InverseDVRmixFE::apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const {
    //HACK apply continuity just to make sure - may be removed eventually
    const_cast<Coefficients*>(&Vec)->makeContinuous();
    inv0->apply(A,Vec,B,Y);
}
void InverseDVRmixFE::applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    Y.makeContinuous();
}

void InverseDVRmixFE::apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{
    const_cast<CoefficientsLocal*>(&Vec)->makeContinuous();
    if(MPIwrapper::Size()>1)DEVABORT("local application not implemented yet");
    inv0->apply(A,Vec,B,Y);
}
void InverseDVRmixFE::applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{
    Y.makeContinuous();
}


