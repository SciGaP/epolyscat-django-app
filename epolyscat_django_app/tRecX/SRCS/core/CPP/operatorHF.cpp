// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorHF.h"
#include "index.h"
#include "basisHF.h"
#include "basisChannel.h"
#include "operatorDefinition.h"

std::shared_ptr<OperatorHartree> OperatorFloorHF::_hfOp;

OperatorFloorHF::OperatorFloorHF(std::string Pot, const Index *IIndex, const Index *JIndex, std::complex<double> Multiplier):
    OperatorFloorNonLin("HF["+Pot+"]")
{
    if(Multiplier!=1.)DEVABORT("for now, cannot have multiplier");
    if(not _hfOp)_hfOp.reset(new OperatorHartree(Pot,IIndex,JIndex,1.));
}

void OperatorFloorHF::updateNonLin(double Time, const Coefficients *C){

}

OperatorHF::OperatorHF(std::string Pot, const Index *IIndex, const Index *JIndex)
    :OperatorNonLin("HF",Pot,IIndex,JIndex)
{
    definition=OperatorDefinition("HartreeFock["+Pot+"]");
    _iChan.reset(new Index());
    _iChan->setAxisName("Channel");
    _iChan->setBasis(new BasisChannel("ChannelHF",new BasisHF(1,IIndex)));
    _iChan->childAdd(new Index(*IIndex));
    _iChan->sizeCompute();
    _applyY.reset(new Coefficients(_iChan->child(0)));
    _applyV.reset(new Coefficients(_iChan->child(0)));

    // definition - somewhat hacky
    std::string coor=_iChan->coordinates();
    std::string def;
    for(int k=0;k<std::count(coor.begin(),coor.end(),'.');k++)def+="<allOnes>";
    _def=OperatorDefinition(def+"<Hartree>");
}

void OperatorHF::update(double Time, const Coefficients *CurrentVec){

    // reset the channel orbital(s)
    dynamic_cast<const BasisHF*>(dynamic_cast<const BasisChannel*>(_iChan->basis())->orbs())->reset({CurrentVec});

    // recompute the hartree operator
    _hf.reset(new OperatorTree("HF",_def,_iChan.get(),_iChan.get()));
}

void OperatorHF::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{

    *_applyV=Vec;
    _hf->child(0)->apply(A,*_applyV,0.,*_applyY);
    Y.axpy(1.,*_applyY,B);
}
