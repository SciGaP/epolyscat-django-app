// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "inverse.h"
#include "useMatrix.h"
#include "tRecXchecks.h"
#include "index.h"

#include "overlapDVR.h"
#include "inverseDvr.h"
#include "inverseFem.h"
#include "inverseDVRmixFE.h"
#include "inverseHybrid.h"
#include "inverseFloors.h"
#include "printOutput.h"

#include "eigenTools.h"

using namespace std;

void Inverse::verify(const OperatorAbstract *Ovr) const{
    if(tRecX::off("inverseOverlap"))return;
    Coefficients one(iIndex);
    Coefficients sOne(iIndex);
    Coefficients res(iIndex);

    UseMatrix Sinv(one.size(),one.size());
    UseMatrix S(one.size(),one.size());
    UseMatrix S_Sinv(one.size(),one.size());
    UseMatrix Sinv_S(one.size(),one.size());
    for(int k=0;k<one.size();k++){
        one.setToZero();
        one.storageData()[k]=1.;
        one.makeContinuous();

        Ovr->apply(1.,one,0.,sOne);
        sOne.makeContinuous();
        S.col(k)=UseMatrix::UseMap(sOne.storageData(),one.size(),1);

        apply(1.,sOne,0.,res);
        Sinv_S.col(k)=UseMatrix::UseMap(res.storageData(),one.size(),1);

        apply(1.,one,0.,sOne);
        Sinv.col(k)=UseMatrix::UseMap(sOne.storageData(),one.size(),1);

        Ovr->apply(1.,sOne,0.,res);
        res.makeContinuous();
        S_Sinv.col(k)=UseMatrix::UseMap(res.storageData(),one.size(),1);
    }
    Sinv_S.purge();
    S_Sinv.purge();
    S.purge();
    S.print("S",0);
    Sinv.purge();
    (S=(S-S.transpose())).purge().print("S",0);
    (Sinv=(Sinv-Sinv.transpose())).purge().print("Sinv",0);
    Sinv_S.print("S_invS",0);
    S_Sinv.print("S_Sinv",0);

}

const Inverse* Inverse::factory(const Index* Idx){
    if(IndexOverlap::inverseOverlap(Idx))return Idx->inverseOverlap(); // already set up
    if(Idx->overlap()==0){
        IndexOverlap::set(Idx);
        PrintOutput::DEVmessage("setup overlap in Inverse::factory");
    }
    Index* idx=const_cast<Index*>(Idx);
    if(idx->hasFloor()){
        idx->setInverseOverlap(new InverseFloors(idx->localOverlap()));
    }
    else if(idx->overlap()->isDiagonal()){
        idx->setInverseOverlap(new InverseDVR(Idx));
    }
    else if(idx->isHybrid()){
        idx->setInverseOverlap(new InverseHybrid(Idx));
    }
    else if(idx->localOverlap()->isBlockDiagonal()) {
        if(idx->localOverlap()==0)DEVABORT("setup local overlap before use Inverse::factory");
        idx->setOverlap(idx->localOverlap());
        idx->setInverseOverlap(new InverseDVRmixFE(idx));
        if(idx->size()<2000){
            // check
            UseMatrix ovr,inv;
            idx->overlap()->matrixAdd(1.,ovr);
            idx->inverseOverlap()->matrixAdd(1.,inv);
            Eigen::Map<Eigen::MatrixXcd>(inv.data(),inv.rows(),inv.cols())*=Eigen::Map<Eigen::MatrixXcd>(ovr.data(),ovr.rows(),ovr.cols());
            if(not inv.isIdentity(1.e-10)){
                Sstr+ovr.str("ovrDirect",1)+Sendl;
                Sstr+ovr.inverse().str("Eigen inverse",1)+Sendl;
                Sstr+inv.str("invDirect",1)+Sendl;
                DEVABORT("inverse FAILED (directly verified): "+idx->overlap()->name);
            }
        }
        else {
            PrintOutput::DEVmessage("inverse OK (directly verified): "+idx->overlap()->name);
        }
    }
    else {
        DEVABORT(Sstr+"no overlap for "+idx+idx->hierarchy());
    }
    return idx->inverseOverlap();
}
void Inverse::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    apply0(A,Vec,B,Y);applyCorrection(A,Vec,B,Y);
}
void Inverse::apply(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{
    apply0(A,Vec,B,Y); applyCorrection(A,Vec,B,Y);
}


