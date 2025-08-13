// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisEigen.h"
#include "operatorTree.h"
#include "operatorDefinition.h"
#include "index.h"
#include "indexNew.h"
#include "eigenSolver.h"
#include "inverse.h"
#include "printOutput.h"
#include "mpiWrapper.h"
#include "parallelOperator.h"
#include "vectorValuedFunction.h"

BasisEigen::BasisEigen(const std::string Def)
    :BasisOrbital("Eigen"){
    std::vector<std::string> part=tools::splitString(Def,':');
    if(part.size()<4)ABORT("need format Eigen:OperDef:RefSubset:Nlow:Size or ...:RefSubset:{0,1,3,..}, got: "+Def);
    _operDef=part[1];
    _refName=part[2];
    if(part.size()==5){
        for(int k=tools::string_to_int(part[3]);k<tools::string_to_int(part[4]);k++)
            _select.push_back(k);
    }
    else
        _select=tools::string_to_intVec("{"+part[3]+"}");
    _orb.resize(_select.size());
    if(_orb.size()==0)DEVABORT("no orbitals set up for BasisEigen: "+Def);
}

BasisEigen* BasisEigen::factory(std::string Def){
    if(Def.find("Eigen:")==0)
        return new BasisEigen(Def);
    if(Def.find("Eigenbasis[")==0){
        return new BasisEigen("Eigen:"+tools::stringInBetween(Def,"[","]"));
    }

    return nullptr;
}

BasisEigen::BasisEigen(const BasisSetDef &Def)
    :BasisEigen("Eigen:"+tools::stringInBetween(Def.funcs,"[","]")+":"+tools::str(Def.lowBound())+":"+tools::str(int(Def.lowBound())+Def.order)){
}

BasisEigen::BasisEigen(std::string OperRef, int Size, int NLow)
    :BasisEigen("Eigen:"+OperRef+":"+tools::str(NLow)+":"+tools::str(NLow+Size)){}

std::string BasisEigen::strDefinition() const {
    return "Eigen:"+_operDef+":"+_refName+":"+tools::str(_select,3,",");
}

void BasisEigen::generateOrbitals(const Index *Idx){

    const Index* idx(Idx);
    if(idx==0)
        idx=BasisOrbital::referenceIndex[_refName];
    if(idx==0){
        idx=BasisOrbital::referenceIndex["main"];
        while(idx!=0 and idx->coordinates()!=_refName)idx=idx->nodeNext();
    }
    if(idx==0)DEVABORT("no index for reference "+_refName+" in hierarchy "+BasisOrbital::referenceIndex["main"]->hierarchy()+
                       "\nset BasisOrbital::referenceIndex[...] before constructing BasisEigen");

    if(not idx->overlap())dynamic_cast<IndexNew*>(const_cast<Index*>(idx))->buildOverlap();

    OperatorTree op("forBasisEigen",OperatorDefinition(_operDef,idx->hierarchy()).str(),idx,idx);
    OperatorTree ov("ovrBasisEigen",OperatorDefinition("<<Overlap>>",idx->hierarchy()).str(),idx,idx);
    ParallelOperator::bcast(&op);
    ParallelOperator::bcast(&ov);
    int kbeg=0,kend=0;
    for(auto k: _select){
        kbeg=std::min(k,kbeg);
        kend=std::max(k+1,kend);
    }
    EigenSolver slv(-DBL_MAX,DBL_MAX,kend,true,true,false,"Lapack");
    slv.withSelect("SmallReal["+tools::str(kend)+"]");

    // things do not work with sparse vectors
    slv.fullVectors().compute(&op,&ov);
    slv.select("SmallReal["+tools::str(kend)+"]");
    slv.orthonormalize();

    _eigenValues.clear();
    for(size_t k=0;k<size();k++){
        _eigenValues.push_back(slv.eigenvalues()[_select[k]]);
        _orb[k]=*slv.rightVectors()[_select[k]];
    }
    slv.verify(1.e-8);
}
