// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisMatOperator.h"

#include "operatorTree.h"
#include "basisOrbital.h"
#include "basisChannel.h"
#include "operatorDefinition.h"
#include "index.h"

BasisMatOperator::BasisMatOperator(std::string Kind)
{
    _operDef=tools::stringInBetween(Kind,"[","]");
    if(_operDef=="")ABORT("must specify operator string for Matrix: operator");
}

void BasisMatOperator::setup(const BasisAbstract *IBas, const BasisAbstract *JBas){

    const BasisOrbital* iBas=dynamic_cast<const BasisOrbital*>(IBas);
    const BasisOrbital* jBas=dynamic_cast<const BasisOrbital*>(JBas);

    if(iBas==0 and jBas==0){
        const BasisChannel* iC=dynamic_cast<const BasisChannel*>(IBas);
        const BasisChannel* jC=dynamic_cast<const BasisChannel*>(JBas);
        if(iC)iBas=iC->orbs();
        if(jC)jBas=jC->orbs();
    }

    if(not iBas)ABORT("Matrix: operator only for Orbital or ChannelHF basis, found lhs basis: "+IBas->str());
    if(not jBas)ABORT("Matrix: operator only for Orbital or ChannelHF basis, found rhs basis: "+JBas->str());
    OperatorTree op("matrix",_operDef,iBas->orbital(0)->idx(),jBas->orbital(0)->idx());
    Coefficients opC(op.iIndex);
    _mat.resize(iBas->size(),jBas->size());
    for(int j=0;j<jBas->size();j++){
        op.apply(1.,*jBas->orbital(j),0.,opC);
        for(int i=0;i<iBas->size();i++){
            _mat(i,j)=iBas->orbital(i)->innerProductUnscaled(&opC);
        }
    }
    EigenTools::purge(_mat);
}
