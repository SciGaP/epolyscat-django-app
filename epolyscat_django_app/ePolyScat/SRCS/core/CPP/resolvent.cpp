// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "resolvent.h"

#include <vector>

#include "operatorTree.h"
#include "str.h"
#include "triFactorLU.h"
#include "index.h"
#include "coefficients.h"
#include "timer.h"
#include "tRecXchecks.h"
#include "blockView.h"
#include "printOutput.h"

TIMER(matrix,resolvent)
TIMER(overlap,resolvent)
TIMER(factorize,resolvent)
TIMER(subtract,resolvent)

Resolvent::~Resolvent(){delete factor;}

Resolvent::Resolvent(const OperatorTree *Op, std::complex<double> Z)
    :OperatorAbstract(Str("1/[","")+Op->name+"-"+Z+"]",Op->iIndex,Op->jIndex),_z(Z),factor(0)
{
    if(Op->isBlockDiagonal()){
        if(Op->descend()->iIndex==Op->iIndex){
            if(Op->childSize()>1)ABORT("DEVELOPER: for now, cannot have composed operators");
            Op=Op->descend();
        }
        for(int k=0;k<Op->childSize();k++){
            childAdd(new Resolvent(Op->child(k),Z));
        }
    }
    else {
        BlockView vOp(Op);

        //note: we should adjust BlockView to handle this automatically
        const OperatorTree* tOv=dynamic_cast<const OperatorTree*>(Op->iIndex->overlap());
        BlockView vOv=tOv?BlockView(tOv):BlockView(Op->iIndex->overlap());

        Eigen::SparseMatrix<std::complex<double> > mOp,mOv;
        vOp.sparseMatrix(mOp,true);
        vOv.sparseMatrix(mOv,true);
        mOp-=mOv*Z;
        solver.analyzePattern(mOp);
        solver.factorize(mOp);
    }
    PrintOutput::message(Sstr+"constructed resolvent for"+Op->name+"at z="+Z);
}


void Resolvent::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{

    if(isLeaf()){
        // contract vector
        Coefficients vec(Vec.idx());
        vec.treeOrderStorage();
        vec=Vec;

        std::vector<unsigned int> idx=iIndex->contractedNumbering();
        UseMatrix vecContract(iIndex->globalLength(),1);

        vec.makeContinuous(sqrt(2.)); // these factors drive me crazy

        for(int k=0;k<vec.size();k++){
            if(idx[k]>vecContract.size())ABORT("bad");
            vecContract.data()[idx[k]]=vec.data()[k];
        }

        Eigen::Map<Eigen::VectorXcd>(vecContract.data(),vecContract.size())=solver.solve(Eigen::Map<Eigen::VectorXcd>(vecContract.data(),vecContract.size()));

        // expand into global storage
        for(int k=0;k<vec.size();k++)vec.data()[k]=vecContract.data()[idx[k]];
        vec.makeContinuous(sqrt(0.5)); // these factors drive me crazy

        Y.scale(B);
        Y.axpy(A,&vec);
    }
    for(int k=0;k<childSize();k++)child(k)->apply(A,*Vec.child(k),B,*Y.child(k));

}

bool Resolvent::verify(const OperatorTree *Op,double Epsilon) const{
    if(tRecX::off("resolvent"))return true;
    Coefficients c(Op->iIndex),d(Op->iIndex);
    c.setToRandom();
    c.makeContinuous();

    apply(1.,c,0.,d);
    Op->apply(-1.,d,1.,c);
    iIndex->overlap()->apply(z(),d,1.,c);
    c.makeContinuous();

    return c.maxCoeff()<Epsilon;
}
