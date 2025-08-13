// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "inverseDvr.h"

#include "mpiWrapper.h"
#include "readInput.h"
#include "operatorTree.h"
#include "coefficients.h"
#include "coefficientsLocal.h"
#include "index.h"
#include "basisIntegrable.h"
#include "str.h"
#include "timer.h"
#include "printOutput.h"
using namespace std;

InverseDVR::InverseDVR(const Index *Idx)
    :Inverse("Inv("+Idx->hierarchy()+")",Idx,Idx),diagonalLocal(0)
{
    // extract diagonal into Coefficients structure, force check
    Coefficients diag(Idx);
    diag.reset(iIndex);
    diag=Idx->localOverlap()->diagonal(true);

    // make continuous
    diag.makeContinuous();

    // get inverse
    _diagonal=new Coefficients(diag.idx(),1.);
    const_cast<Coefficients*>(_diagonal)->treeOrderStorage();
    const_cast<Coefficients*>(_diagonal)->cwiseDivide(diag);

    // create hierachy until first continuity
    if(Idx->continuity()==Index::npos)
        for(size_t k=0;k<_diagonal->childSize();k++)
            childAdd(new InverseDVR(_diagonal->child(k)));
    if(!Idx->parent() and _diagonal->isNan()){
        DEVABORT("_diagonal.isNan");
    }
}

InverseDVR::InverseDVR(const Coefficients *Diagonal)
    :Inverse("Inv("+Diagonal->idx()->hierarchy()+")",Diagonal->idx(),Diagonal->idx()),
      _diagonal(Diagonal),diagonalLocal(0)
{
    // descend hierarchy until first continuity level
    if(Diagonal->parent()==0 or Diagonal->parent()->idx()->continuity()==Index::npos)
        for(size_t k=0;k<Diagonal->childSize();k++)
            Diagonal->idx()->child(k)->setInverseOverlap(new InverseDVR(Diagonal->child(k)));
    if(!Diagonal->parent() and _diagonal->isNan())DEVABORT("diag");
}

void InverseDVR::apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    if(A!=1.)ABORT("A!=1 not implemented");
    if(B!=0.)*tempLHS()=Y;
    Y=Vec;
    Y.cwiseMultiply(*_diagonal);
    if(B!=0.)Y.axpy(B,tempLHS());
}

void InverseDVR::parallelSetup() const{
    const_cast<InverseDVR*>(this)->diagonalLocal=CoefficientsLocal::view(const_cast<Coefficients*>(_diagonal));
}

void InverseDVR::apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{
    if(A!=1.)ABORT("A!=1 not implemented");
    if(B!=0.)ABORT("B!=0 not implemented");
    if(diagonalLocal==0)ABORT("call parallelSetup() before use with CoefficientsLocal");
    Y=Vec;
    Y.cwiseMultiply(*diagonalLocal);
}
