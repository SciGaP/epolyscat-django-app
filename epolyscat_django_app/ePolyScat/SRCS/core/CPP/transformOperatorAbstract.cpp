// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "transformOperatorAbstract.h"

#include "operatorFloor.h"
#include "coefficients.h"
#include "index.h"
#include "operatorTree.h"
#include "operatorAbstract.h"
#include "operatorAbstractProduct.h"


TransformOperatorAbstract::TransformOperatorAbstract(const OperatorAbstract* Op, const OperatorAbstract* LeftTrafo, const OperatorAbstract* RightTrafo):
    op(Op), leftTrafo(LeftTrafo), rightTrafo(RightTrafo), transformed(0){
    
    if(op->iIndex!=rightTrafo->jIndex or op->jIndex!=leftTrafo->iIndex) ABORT("Index mismatch!");
}


void TransformOperatorAbstract::build(OperatorTree* optree, OperatorAbstract* op, UseMatrix* mat, const Index* matIIndex, const Index* matJIndex){
    UseMatrix tmp;

    if(mat==0 and (optree->iIndex->sizeStored()*optree->jIndex->sizeStored()<4000*4000 or (optree->iIndex->hasFloor() and optree->jIndex->hasFloor()))){

        std::vector<std::complex<double> > result_;

        op->subMatrix(
                result_,
                optree->iIndex, 
                optree->jIndex        
        );
        
        tmp=UseMatrix(UseMatrix::UseMap(result_.data(), optree->iIndex->sizeStored(), optree->jIndex->sizeStored()));
        tmp.purge();
        
        mat=&tmp;
        matIIndex=optree->iIndex;
        matJIndex=optree->jIndex;
    }

    if(mat!=0 and mat->isZero()) return;

    if(optree->iIndex->hasFloor() and optree->jIndex->hasFloor()){
        UseMatrix result = mat->block(
            optree->iIndex->posIndex(matIIndex),
            optree->jIndex->posIndex(matJIndex),
            optree->iIndex->sizeStored(),
            optree->jIndex->sizeStored()
        );
        optree->oFloor = OperatorFloor::factory(std::vector<const UseMatrix*>(1, &result), optree->name);

    }else if(optree->iIndex->hasFloor() and not optree->jIndex->hasFloor()){
        for(int i=0; i<optree->jIndex->childSize(); i++){
            OperatorTree* c = new OperatorTree(optree->name, optree->iIndex, optree->jIndex->child(i));
            build(c, op, mat, matIIndex, matJIndex);
            if(c->isZero()){
                delete c;
            }else{
                optree->childAdd(c);
            }
        }
    }else if(not optree->iIndex->hasFloor() and optree->jIndex->hasFloor()){
        for(int i=0; i<optree->iIndex->childSize(); i++){
            OperatorTree* c = new OperatorTree(optree->name, optree->iIndex->child(i), optree->jIndex);
            build(c, op, mat, matIIndex, matJIndex);
            if(c->isZero()){
                delete c;
            }else{
                optree->childAdd(c);
            }
        }
    }else{
        int l=optree->iIndex->depth();
        bool diagonal=false;
        for(int i=0; i<optimize__diagonalLevels.size(); i++) if(optimize__diagonalLevels[i]==l) diagonal=true;
        
        for(int i=0; i<optree->iIndex->childSize(); i++){
            for(int j=0; j<optree->jIndex->childSize(); j++){
                if(diagonal and i!=j) continue;

                OperatorTree* c = new OperatorTree(optree->name, optree->iIndex->child(i), optree->jIndex->child(j));
                build(c, op, mat, matIIndex, matJIndex);
                if(c->isZero()){
                    delete c;
                }else{
                    optree->childAdd(c);
                }
            }
        }

    }
}

void TransformOperatorAbstract::transform(){
    transformed = new OperatorTree(op->name, rightTrafo->iIndex, leftTrafo->jIndex);

    std::vector<const OperatorAbstract*> ops;
    ops.push_back(rightTrafo);
    ops.push_back(op);
    ops.push_back(leftTrafo);
    OperatorAbstractProduct product("Transformed("+op->name+")",ops);

    std::cerr<<"Building transformed optree... ";
    build(transformed, &product);
    std::cerr<<"done"<<std::endl;
}
    
OperatorAbstract* TransformOperatorAbstract::getTransformedWithTransformations(){

    std::vector<const OperatorAbstract*> ops;
    ops.push_back(leftTrafo);
    ops.push_back(transformed);
    ops.push_back(rightTrafo);
    return new OperatorAbstractProduct(transformed->name,ops);
}


bool TransformOperatorAbstract::makeDiagonal(std::map<std::string, double>& report, double eps){
    if(transformed==0) transform();
    
    // OPTIONAL
    report["transformedDim"]=rightTrafo->iIndex->sizeStored();
    report["originalDim"]=rightTrafo->jIndex->sizeStored();

    // Create the projector
    std::vector<const OperatorAbstract*> P_tmp0;
    P_tmp0.push_back(rightTrafo);
    P_tmp0.push_back(leftTrafo);
    OperatorAbstractProduct P_tmp1("P", P_tmp0);

    // Get the projectors matrix P
    UseMatrix P_tmp2;
    P_tmp1.matrix(P_tmp2);
    UseMatrix P(P_tmp2);

    // Create matrix A
    UseMatrix PPdagger = P*P.adjoint();
    UseMatrix PdaggerP = P.adjoint()*P;
    UseMatrix A(P.rows(), P.cols());

    for(int i=0; i<A.rows(); i++){
        for(int j=0; j<A.cols(); j++){
            A(i,j) = PdaggerP(i,j)*PPdagger(j,i);
        }
    }

    // Diagonalise A and create the pseudoinverse matrix Ainv
    UseMatrix Aeigenval, Aeigenvec, dummy;
    A.eigen(Aeigenval, Aeigenvec);
    A.eigenOrthonormalize(Aeigenval, Aeigenvec, dummy);
    UseMatrix AinvDiagonal=UseMatrix::Zero(Aeigenval.rows(), Aeigenval.rows());
    int AZeroEigenvals=0;
    for(int i=0; i<Aeigenval.rows(); i++){
        if(std::abs(Aeigenval(i).complex())>eps){
            AinvDiagonal(i,i)=1./Aeigenval(i);
        }else AZeroEigenvals++;
    }
    report["AZeroEigenvals"]=AZeroEigenvals;

    UseMatrix Ainv = Aeigenvec * AinvDiagonal * Aeigenvec.adjoint();

    // Create vector b
    UseMatrix V_tmp0;
    transformed->matrix(V_tmp0);
    UseMatrix V(V_tmp0);

    UseMatrix PdaggerVPdagger = P.adjoint() * V * P.adjoint();
    UseMatrix b(PdaggerVPdagger.rows(), 1);
    for(int i=0; i<b.rows(); i++){
        b(i)=PdaggerVPdagger(i,i);
    }

    // OPTIONAL
    // Investigate b\in\mathcal A
    UseMatrix tildeB = Aeigenvec.adjoint() * b;
    int bInAZeroElements=0;
    for(int i=0; i<tildeB.rows(); i++){
        if(std::abs(Aeigenval(i).complex())<eps and std::abs(tildeB(i).complex())>eps) bInAZeroElements++;
    }
    report["bInAZeroElements"]=bInAZeroElements;

    // OPTIONAL
    // Calculate tildeC
    double normsqV=0;
    for(int i=0; i<V.rows(); i++){
        for(int j=0; j<V.cols(); j++){
            normsqV+=std::pow(std::abs(V(i,j).complex()),2);
        }
    }

    UseMatrix tildeC_tmp0 = b.adjoint() * Ainv * b;
    report["tildeC"] = tildeC_tmp0(0,0).real() - normsqV;

    // Calculate the new diagonal matrix
    UseMatrix diagonal = Ainv*b;
    UseMatrix diagonalMat=UseMatrix::Zero(diagonal.rows(), diagonal.rows());
    for(int i=0; i<diagonal.rows(); i++){
        diagonalMat(i,i)=diagonal(i);
    }

    OperatorFloor* floor = OperatorFloor::factory(std::vector<const UseMatrix*>(1, &diagonalMat), transformed->name);
    OperatorTree* transformedDiagonal = new OperatorTree(transformed->name, "",  transformed->iIndex, transformed->jIndex, floor);


    // Check if \norm{leftTrafo*transfromedDiagonal*rightTrafo - op}<eps
    UseMatrix ownMat;
    UseMatrix opMat;

    std::vector<const OperatorAbstract*> own_tmp0;
    own_tmp0.push_back(leftTrafo);
    own_tmp0.push_back(transformedDiagonal);
    own_tmp0.push_back(rightTrafo);
    OperatorAbstractProduct own_tmp1("", own_tmp0);
    own_tmp1.matrix(ownMat);

    op->matrix(opMat);

    double normsqDiff=0.;
    double normsqOp=0.;
    for(int i=0; i<opMat.rows(); i++){
        for(int j=0; j<opMat.cols(); j++){
            normsqDiff+=std::pow(std::abs(ownMat(i,j).complex()-opMat(i,j).complex()),2);
            normsqOp+=std::pow(std::abs(opMat(i,j).complex()),2);
        }
    }
    report["diff"] = std::sqrt(normsqDiff/normsqOp);

    if(std::sqrt(normsqDiff/normsqOp)<eps){
        delete transformed;
        transformed = transformedDiagonal;
        return true;
    }else{
        delete transformedDiagonal;
        return false;
    }
}


void TransformOperatorAbstract::check(){
    Coefficients c(op->jIndex);
    Coefficients cRight(transformed->jIndex);
    Coefficients cLeft(transformed->iIndex);
    Coefficients c1(op->iIndex);
    Coefficients c2(op->iIndex);

    c.treeOrderStorage();
    cRight.treeOrderStorage();
    cLeft.treeOrderStorage();
    c1.treeOrderStorage();
    c2.treeOrderStorage();

    for(int i=0; i<100; i++){
        c.setToRandom();
        op->apply(1., c, 0., c1);

        rightTrafo->apply(1., c, 0., cRight);
        transformed->apply(1., cRight, 0., cLeft);
        leftTrafo->apply(1., cLeft, 0., c2);
        
        c1-=c2;
        if(c1.maxCoeff()>1.e-9) std::cerr<<("WARNING! Transfom din't work: error="+std::to_string(c1.maxCoeff())+"\n");
    }
}











