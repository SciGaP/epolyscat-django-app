// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorTucker.h"
#include "index.h"
#include "coefficients.h"
#include "transformOperatorAbstract.h"

#include <vector>
#include <complex>


void OperatorTucker::TensorProduct::fillIndicesRec(int* fillIndices, int offset, std::vector<int> dims, std::vector<int> filledDims){
    int dim=dims[0];
    if(dims.size()==1){
        for(int i=0; i<dim; i++){
            fillIndices[i]=offset+i;
        }
    }else{
        dims.erase(dims.begin());
        filledDims.erase(filledDims.begin());

        int scale=1;
        int o=1;
        for(int j=0; j<filledDims.size(); j++) scale*=filledDims[j];
        for(int j=0; j<dims.size(); j++) o*=dims[j];
        for(int i=0; i<dim; i++){
            fillIndicesRec(fillIndices+i*o, offset+i*scale, dims, filledDims);
        }

    }
}
    

OperatorTucker::TensorProduct::TensorProduct(std::string Name, const Index* IIndex, const Index* JIndex, std::vector<UseMatrix*> Mats):
    OperatorAbstract(Name, IIndex, JIndex), mats(Mats){

    std::vector<int> iDims;
    std::vector<int> jDims;
    std::vector<int> filledDims;

    int tempSize = 1;
    for(int i=0; i<mats.size(); i++){
        iDims.push_back(mats[i]->rows());
        jDims.push_back(mats[i]->cols());
        filledDims.push_back(std::max(mats[i]->rows(), mats[i]->cols()));
        tempSize*=std::max(mats[i]->rows(), mats[i]->cols());
    }

    temporary1 = new std::vector<std::complex<double> >(tempSize);
    temporary2 = new std::vector<std::complex<double> >(tempSize);

    fillIndicesI=std::vector<int>(iIndex->sizeStored());
    fillIndicesJ=std::vector<int>(jIndex->sizeStored());

    // TODO! Probably this doesn't work, if iDims!=jDims; never tested, at least
    fillIndicesRec(fillIndicesI.data(), 0, iDims, filledDims);
    fillIndicesRec(fillIndicesJ.data(), 0, jDims, filledDims);

}



void OperatorTucker::TensorProduct::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    TensorProduct* This = const_cast<TensorProduct*>(this);

    // Feed temporary1<-Fill[Vec]
    for(int j=0; j<temporary1->size(); j++) (*temporary1)[j]=0.;
    std::vector<std::complex<double>* > pVec;
    const_cast<Coefficients&>(Vec).pointerToC(pVec);
    for(int j=0; j<fillIndicesJ.size(); j++) (*temporary1)[fillIndicesJ[j]]=*pVec[j];

    // Apply temporary1<-(FillSquare[mats[0]] \tensor Identity \tensor ...)(Identity \tensor FillSquare[mats[1]] \tensor Identity ...) ... temporary1
    for(int i=0; i<mats.size(); i++){
        int identityDimBeforeFactor=1;
        int identityDimAfterFactor=1;
        for(int j=0; j<i; j++) identityDimBeforeFactor*=std::max(mats[j]->rows(), mats[j]->cols());
        for(int j=i+1; j<mats.size(); j++) identityDimAfterFactor*=std::max(mats[j]->rows(), mats[j]->cols());


        for(int j=0; j<temporary2->size(); j++) (*This->temporary2)[j]=0.;
        
        for(int b=0; b<identityDimBeforeFactor; b++){
            for(int m=0; m<mats[i]->rows(); m++){
                for(int n=0; n<mats[i]->cols(); n++){
                    for(int a=0; a<identityDimAfterFactor; a++){
                        (*This->temporary2)[b*identityDimAfterFactor*std::max(mats[i]->rows(), mats[i]->cols()) + m*identityDimAfterFactor + a]+=(*mats[i])(m,n)*(*This->temporary1)[b*identityDimAfterFactor*std::max(mats[i]->rows(), mats[i]->cols()) + n*identityDimAfterFactor + a];
                    }
                }
            }
        }

        std::vector<std::complex<double> >* tmp = temporary2;
        This->temporary2 = temporary1;
        This->temporary1 = tmp;
    }


    Y*=B;

    // Feed Y<-Y + A*Unfill[temporary1] 
    std::vector<std::complex<double>* > pY;
    Y.pointerToC(pY);
    for(int j=0; j<fillIndicesI.size(); j++) *pY[j]+=A*(*temporary1)[fillIndicesI[j]];
}

long OperatorTucker::TensorProduct::applyCount() const{
    long result=0;
    for(int i=0; i<mats.size(); i++){
        int identityDimBeforeFactor=1;
        int identityDimAfterFactor=1;
        for(int j=0; j<i; j++) identityDimBeforeFactor*=std::max(mats[j]->rows(), mats[j]->cols());
        for(int j=i+1; j<mats.size(); j++) identityDimAfterFactor*=std::max(mats[j]->rows(), mats[j]->cols());
        
        result+=identityDimBeforeFactor*identityDimAfterFactor*mats[i]->rows()*mats[i]->cols();
    }
    return result;

}


void OperatorTucker::ensureTensorProduct(const Index* idx, const Index* base, std::vector<int>& dims){
    if(idx->childSize()==0) return;

    int i = idx->depth()-base->depth();

    if(i>=dims.size()){
        dims.push_back(idx->childSize());
    }else{
        if(idx->childSize()!=dims[i]) ABORT("Not tensor product structure");
    }

    for(int i=0; i<idx->childSize(); i++) ensureTensorProduct(idx->child(i), base, dims);

}


void OperatorTucker::calculateDensityMatrix(UseMatrix& mat, const OperatorAbstract* op, bool isLHS, int order){
    
    std::vector<const Index*> parents;
    if(isLHS) op->iIndex->childrenAtDepth(order, parents);
    else      op->jIndex->childrenAtDepth(order, parents);

    mat=UseMatrix(parents[0]->childSize(), parents[0]->childSize());
    for(int m=0; m<mat.rows(); m++){
        for(int n=0; n<mat.cols(); n++){
            mat(m,n)=0.;
            for(int i=0; i<parents.size(); i++){
                const Index* mIdx = parents[i]->child(m);
                const Index* nIdx = parents[i]->child(n);

                std::vector<std::complex<double> > mMat;
                std::vector<std::complex<double> > nMat;

                if(isLHS){
                    op->subMatrix(mMat, mIdx, op->jIndex);
                    op->subMatrix(nMat, nIdx, op->jIndex);
                }else{
                    op->subMatrix(mMat, op->iIndex, mIdx);
                    op->subMatrix(nMat, op->iIndex, nIdx);
                }

                if(mMat.size()!=nMat.size()) ABORT("Sizes mismatch");
                for(int i=0; i<mMat.size(); i++) mat(m,n)+=std::conj(mMat[i])*nMat[i];

            }
        }
    }
}


void OperatorTucker::createIndex(Index* idx, std::vector<int> hierarchy){
    if(hierarchy.size()==0) return;

    for(int i=0; i<hierarchy[0]; i++){
        idx->childAdd(new Index());
    }

    hierarchy.erase(hierarchy.begin());

    for(int i=0; i<idx->childSize(); i++){
        createIndex(idx->child(i), hierarchy);
    }
}


OperatorTucker* OperatorTucker::truncate(const OperatorAbstract* base){
    std::vector<int> iDims;
    ensureTensorProduct(base->iIndex, base->iIndex, iDims);
    std::vector<int> jDims;
    ensureTensorProduct(base->jIndex, base->jIndex, jDims);
    
    std::vector<UseMatrix*> iTransformations;
    std::vector<UseMatrix*> jTransformations;

    for(int i=0; i<iDims.size(); i++){
        std::cerr<<"CALCULATING LHS DENSITY MATRIX "<<i<<std::endl;
        UseMatrix densityMatrix;
        calculateDensityMatrix(densityMatrix, base, true, i);
        

        UseMatrix weights;
        UseMatrix transformation;

        if(densityMatrix.rows()>1 or densityMatrix.cols()>1){
            UseMatrix ovr;

            densityMatrix.eigen(weights, transformation, ovr);
            densityMatrix.eigenOrthonormalize(weights, transformation, ovr);
        }else{
            weights = UseMatrix(1,1);
            weights(0)=1.;
            transformation = UseMatrix(1,1);
            transformation(0,0)=1.;
        }

        weights.print();

        UseMatrix check(transformation);
        check*=check.adjoint();
        if(not check.isIdentity(1.e-9)) ABORT("Density matrix eigen base is not ONB");
        
        iTransformations.push_back(new UseMatrix(transformation));
    }
    
    for(int i=0; i<jDims.size(); i++){
        std::cerr<<"CALCULATING RHS DENSITY MATRIX "<<i<<std::endl;
        UseMatrix densityMatrix;
        calculateDensityMatrix(densityMatrix, base, false, i);
        
        UseMatrix weights;
        UseMatrix transformation;
        
        if(densityMatrix.rows()>1 or densityMatrix.cols()>1){
            UseMatrix ovr;

            densityMatrix.eigen(weights, transformation, ovr);
            densityMatrix.eigenOrthonormalize(weights, transformation, ovr);
        }else{
            weights = UseMatrix(1,1);
            weights(0)=1.;
            transformation = UseMatrix(1,1);
            transformation(0,0)=1.;
        }
        
        weights.print();
       
        UseMatrix transformationAdj(transformation.adjoint());
        
        UseMatrix check(transformationAdj);
        check*=check.adjoint();
        if(not check.isIdentity(1.e-9)) ABORT("Density matrix eigen base is not ONB");
        
        jTransformations.push_back(new UseMatrix(transformationAdj));
    }

    Index* iInnerIndex = new Index();
    Index* jInnerIndex = new Index();

    std::vector<int> iHierarchy;
    std::vector<int> jHierarchy;
    for(int i=0; i<iTransformations.size(); i++) iHierarchy.push_back(iTransformations[i]->cols());
    for(int i=0; i<jTransformations.size(); i++) jHierarchy.push_back(jTransformations[i]->rows());

    createIndex(iInnerIndex, iHierarchy);
    createIndex(jInnerIndex, jHierarchy);

    iInnerIndex->resetFloor(base->iIndex->firstFloor()->depth() - base->iIndex->depth());
    jInnerIndex->resetFloor(base->jIndex->firstFloor()->depth() - base->jIndex->depth());

    iInnerIndex->sizeCompute();
    jInnerIndex->sizeCompute();
    
    TensorProduct* mapToInner = new TensorProduct("MapToInner("+base->name+")", jInnerIndex, base->jIndex, jTransformations);
    TensorProduct* mapFromInner = new TensorProduct("MapFromInner("+base->name+")", base->iIndex, iInnerIndex, iTransformations);

    std::vector<UseMatrix*> iTransformationsAdj;
    std::vector<UseMatrix*> jTransformationsAdj;

    for(int i=0; i<iTransformations.size(); i++) iTransformationsAdj.push_back(new UseMatrix(iTransformations[i]->adjoint()));
    for(int i=0; i<jTransformations.size(); i++) jTransformationsAdj.push_back(new UseMatrix(jTransformations[i]->adjoint()));

    TensorProduct mapToInnerAdj("MapToInnerAdj("+base->name+")", base->jIndex, jInnerIndex, jTransformationsAdj); 
    TensorProduct mapFromInnerAdj("MapFromInnerAdj("+base->name+")", iInnerIndex, base->iIndex, iTransformationsAdj); 

    std::cerr<<"CALCULATING CORE TENSOR"<<std::endl;

    TransformOperatorAbstract transform(base, &mapToInnerAdj, &mapFromInnerAdj);
    transform.transform();
    transform.check();
    OperatorAbstract* core = transform.getTransformed();

    UseMatrix mat;
    core->matrix(mat);
    mat.print();
    mat.show();

    std::cerr<<"MAP FROM INNER: "<<mapFromInner->applyCount()<<std::endl;
    std::cerr<<"MAP TO INNER: "<<mapToInner->applyCount()<<std::endl;
    std::cerr<<"CORE: "<<core->applyCount()<<std::endl;
    std::cerr<<"BASE: "<<base->applyCount()<<std::endl;

    std::vector<const OperatorAbstract*> maps;
    maps.push_back(mapFromInner);
    maps.push_back(core);
    maps.push_back(mapToInner);
    OperatorTucker* result = new OperatorTucker("TuckerRep("+base->name+")", maps);
    result->check(base);
    return result;
}


OperatorTucker::OperatorTucker(std::string name, std::vector<const OperatorAbstract*> maps):
    OperatorAbstractProduct(name, maps){

    mapFromInner = maps[0];
    core = maps[1];
    mapToInner = maps[2];
}

void OperatorTucker::check(const OperatorAbstract* base){
    Coefficients c(base->jIndex);
    Coefficients c1(base->iIndex);

    for(int i=0; i<100; i++){
        c.setToRandom();
        apply(1.,c,0.,c1);
        base->apply(1.,c,-1.,c1);
        if(not c1.isZero(1.e-9)) ABORT("Truncation did not work");
    }
}













