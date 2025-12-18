// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "inverseHybrid.h"

#include "coefficientsGlobal.h"
#include "discretizationHybrid.h"
#include "printOutput.h"
#include "operatorVectors.h"
#include "mpiWrapper.h"
#include "parallel.h"

static std::unique_ptr<CoefficientsGlobal> gY, gX;

InverseHybrid::InverseHybrid(const Index* Idx)
    :Inverse("Inv("+Idx->axisName()+")",Idx,Idx),sAB(0){

    if(not Idx->isHybrid() or Idx->childSize()!=2)
        DEVABORT("for now only two hydbrid components on top");

    if(Idx->child(0)->size()>Idx->child(1)->size())
        PrintOutput::warning("sort summands such that first is smaller");

    // inverses for sub-blocks may not exist yet - set recursively
    for(size_t k=0;k<Idx->childSize();k++){
        if(!Idx->child(k)->inverseOverlap())
            Idx->child(k)->setInverseOverlap(Inverse::factory(Idx->child(k)));
    }

    sAinv=Idx->child(0)->inverseOverlap();
    sBinv=Idx->child(1)->inverseOverlap();

    // overlap block, local pointers
    const OperatorTree*block01=0,*block10=0,*block00=0;
    for(const OperatorTree* o=dynamic_cast<const OperatorTree*>(Idx->overlap());o!=0;o=o->nodeNext()){
        if(o->iIndex==Idx->child(0) and o->jIndex==Idx->child(0))block00=o;
        if(o->iIndex==Idx->child(0) and o->jIndex==Idx->child(1))block01=o;
        if(o->iIndex==Idx->child(1) and o->jIndex==Idx->child(0))block10=o;
        if(block00 and block01 and block10)break;
    }
    sAB=block01;

    aVec.reset(new Coefficients(Idx->child(0)));
    bVec.reset(new Coefficients(Idx->child(1)));

    zInv=Eigen::MatrixXcd::Zero(aVec->size(),aVec->size());
    if(!block01 and !block10 and Idx->axisSubset()!="Subspace&Complement"){
        PrintOutput::warning("hybrid blocks seem to be orthogonal for "+Idx->axisName());
    }
    Coefficients aJ(aVec->idx());
    Coefficients bTmp(bVec->idx());
    if(sAB!=0)sInvBA.reset(new OperatorVectors("SbInvBA",block10->iIndex,block10->jIndex));
    for(size_t j=0;j<aVec->size();j++){
        aJ.setToZero();
        aJ.data()[j]=1.;// j'th unit vector ej
        block00->apply(1.,aJ,0.,*aVec);  // a = Sa ej
        if(sAB!=0){
            block10->apply(1.,aJ,0.,*bVec); // b = C ej
            sBinv->apply(1.,*bVec,0.,bTmp);
            sInvBA->insertColumn(j,bTmp);
            sAB->apply(-1.,bTmp,1.,*aVec); // a <- a - C^H Sb^-1 C ej = [Sa - C^H Sb^-1 C] ej
        }
        zInv.col(j)=Eigen::Map<Eigen::MatrixXcd>(aVec->data(),zInv.rows(),1);
    }
    zInv=zInv.inverse();
//    Eigen::JacobiSVD<Eigen::MatrixXcd> jDec(zInv,Eigen::ComputeThinU | Eigen::ComputeThinV);
//    zInv.setZero();
//    for(int k=0;k<jDec.singularValues().size();k++){
//        if(std::abs(jDec.singularValues()(k))>1.e-12){
//            zInv=jDec.matrixU().col(k)*jDec.singularValues()(k)*jDec.matrixV().col(k).adjoint();
//        }
//    }
}

void InverseHybrid::parallelSetup() const{
    if(MPIwrapper::Size()>1) {//ABORT("for now, only single processor");
        gX.reset(new CoefficientsGlobal(jIndex));
        gY.reset(new CoefficientsGlobal(iIndex));
    }
}
void InverseHybrid::apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    DEVABORT("not implemented");
}
void InverseHybrid::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    // use Woodbury formula, see tsurff.pdf for the algorithm
    if(Y.childSize()!=2)ABORT("hybrid basis must have exactly 2 components");

    Y.scale(B); // res = B*Y;
    if(A!=1.)ABORT("not for A!=1");

    *aVec=*Vec.child(0); // a = va
    sBinv->apply(1.,*Vec.child(1),0.,*bVec); // b = Sb^-1 vb
    if(sAB!=0)sAB->apply(-1.,*bVec,1.,*aVec); // a = va - C^T b
    Eigen::Map<Eigen::VectorXcd>(aVec->data(),zInv.rows())=zInv*Eigen::Map<Eigen::VectorXcd>(aVec->data(),zInv.rows());// a <- Z^-1 a
    if(sAB!=0)sInvBA->apply(-1.,*aVec,1.,*Y.child(1));  // res=res - Sb^-1 C a = res - Sb^-1 C Z^-1 (va - C^T Sb^-1 vb)

    Y.child(0)->operator+=(*aVec); // res=res + a = res + Z^-1(va-C^T Sb^-1 b)
    Y.child(1)->operator+=(*bVec); // res=res + b = res + Sb^-1 vb
}

void InverseHybrid::apply(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{
    //HACK safe and slow - gather/apply/scatter
    Coefficients* yPtr;
    const Coefficients* vecPtr;
    if(MPIwrapper::Size() > 1) {
        Parallel::allGather(gX.get(), const_cast<CoefficientsLocal*>(dynamic_cast<const CoefficientsLocal*>(&Vec)));
        vecPtr = gX.get();
        yPtr = gY.get();
    }
    else {
        yPtr = &Y;
        vecPtr = &Vec;
    }
    apply(A,*dynamic_cast<const Coefficients*>(vecPtr),B,*dynamic_cast<Coefficients*>(yPtr));
    if(MPIwrapper::Size() > 1)
        Parallel::scatter(gY.get(), dynamic_cast<CoefficientsLocal*>(&Y), MPIwrapper::master());
}
void InverseHybrid::applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    DEVABORT("not implemented");
}
void InverseHybrid::apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{
    DEVABORT("not implemented");
}
void InverseHybrid::applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{
    DEVABORT("not implemented");
}
