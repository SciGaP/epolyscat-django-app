// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorSubspace.h"

#include "projectSubspace.h"
#include "operatorTree.h"
#include "basisOrbital.h"
#include "eigenSolver.h"
#include "inverse.h"

OperatorSubspace::OperatorSubspace(const OperatorAbstract* Op,std::shared_ptr<const ProjectSubspace> Proj)
    :OperatorTree(Op->name,Op->iIndex,Op->jIndex),_proj(Proj),_op(Op){
    if(Proj!=0 and (Proj->iIndex!=Op->iIndex or Proj->iIndex!=Op->jIndex))
        ABORT("Operator and projector are not on the same space");
    definition=OperatorDefinition(Op->def());
    _C.reset(new Coefficients(iIndex));
    _D.reset(new Coefficients(iIndex));
}

OperatorSubspace::OperatorSubspace(const OperatorAbstract *Op):OperatorSubspace(Op,0)
{
    if(iIndex->axisSubset()!="Subspace&Complement")
        ABORT("cannot infer a projection for index:\n"+iIndex->str());

    const BasisOrbital * b=dynamic_cast<const BasisOrbital*>(iIndex->child(0)->basis());
    if(b==0)ABORT("Subspace does not seem to have orbital basis: "+iIndex->child(0)->strNode());

    std::vector<Coefficients*>vecs,dual;
    Coefficients v(iIndex->child(1));
    for(size_t k=0;k<b->size();k++){
        vecs.push_back(new Coefficients(iIndex));
        dual.push_back(new Coefficients(iIndex));
        *vecs.back()->child(1)=*b->orbital(k);
        v=*b->orbital(k);
        iIndex->child(1)->overlap()->apply(1.,v,0.,*dual.back()->child(1));
    }
    _proj.reset(new ProjectSubspace(vecs,dual));
}

void OperatorSubspace::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    if(_proj==0)DEVABORT("incomplete setup - no projector defined");
    *_C=Vec;
    _proj->apply(-1.,Vec,1.,*_C);
    _op->apply(A,*_C,0.,*_D);
    _op->idx()->inverseOverlap()->apply(1.,*_D,0,*_C);
    _proj->apply(-1.,*_D,1.,*_C);
    _op->idx()->overlap()->apply(1.,*_C,B,Y);
}

