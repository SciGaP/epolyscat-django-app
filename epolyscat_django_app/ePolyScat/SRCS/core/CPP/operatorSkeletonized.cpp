// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
//
//  operatorSkeletonized.cpp
//  tRecX
//
//  Created by Jonas Bucher on 19.04.17.
//  Copyright © 2017 youdontneedtoknow. All rights reserved.
//

#include <vector>
#include <complex>

#include "operatorSkeletonized.h"
#include "tools.h"

#include "index.h"
#include "coefficients.h"
#include "abort.h"

#include "str.h"

#include "eigenNames.h"


void OperatorSkeletonized::apply(std::complex<double> alpha, const Coefficients &x, std::complex<double> beta, Coefficients &y) const{
    if(x.idx()!=idx) ABORT("rhs indices don't match");
    if(y.idx()!=idx) ABORT("lhs indices don't match");
    if(S==0) ABORT("OperatorSkeletonized::apply called on non-root");
    
    if(beta==0.)y.setToZero();
    else if(beta!=1.)y*=beta;
    beta=1.;
    
    
    /*
     * Compress x
     */
    std::vector<std::complex<double> > compressed;
    compressed.resize(S->cols());
    
    unsigned int c=0;
    for(unsigned int i=0; i<childSize(); i++){
        child(i)->compressStoringDiagonals(compressed.data()+c, *x.child(i));
        c+=child(i)->R->rows();
    }
    
    /*
     * Apply S to compressed x
     */
    std::vector<std::complex<double> > applied;
    applied.resize(S->rows());
    Map<VectorXcd>(applied.data(), S->rows())=(*S)*Map<VectorXcd>(compressed.data(), S->cols());
    
    /*
     * Take care of alpha
     */
    for(unsigned int i=0; i<applied.size();i++) applied[i]*=alpha;
    
    /*
     * Decompress into y
     */
    c=0;
    for(unsigned int i=0; i<childSize(); i++){
        child(i)->decompressRestoringDiagonals(applied.data()+c, *y.child(i));
        c+=child(i)->L->cols();
    }
    
}

void OperatorSkeletonized::compressStoringDiagonals(std::complex<double>* target, const Coefficients& x) const{
    if(x.idx()!=idx) ABORT("rhs indices don't match");
    
    std::complex<double>* source = temporaryStorage->data();
    
    if(childSize()==0){
        
        std::vector<std::complex<double>* > pC;
        const_cast<Coefficients&>(x).pointerToC(pC);
        for(int i=0; i<R->cols(); i++){
            (*temporaryStorage)[i]=*pC[i];
        }
        
    }else{
        
        unsigned int c=0;
        for(unsigned int i=0; i<childSize(); i++){
            child(i)->compressStoringDiagonals(source+c, *x.child(i));
            c+=child(i)->R->rows();
        }
        
    }
    
    Map<VectorXcd>(target, R->rows())=(*R)*Map<VectorXcd>(source,R->cols());
    Map<VectorXcd>(source, D->rows())=(*D)*Map<VectorXcd>(source,D->cols());
}

void OperatorSkeletonized::decompressRestoringDiagonals(std::complex<double> *source, Coefficients& y) const{
    if(y.idx()!=idx) ABORT("lhs indices don't match");
    
    Map<VectorXcd>(temporaryStorage->data(),L->rows())+=(*L)*Map<VectorXcd>(source,L->cols());
    
    if(childSize()==0){
        
        std::vector<std::complex<double>* > pC;
        y.pointerToC(pC);
        for(int i=0; i<R->cols(); i++){
            *(pC[i])+=(*temporaryStorage)[i];
        }
        
    }else{
        
        unsigned int c=0;
        for(unsigned int i=0; i<childSize(); i++){
            child(i)->decompressRestoringDiagonals(temporaryStorage->data()+c, *y.child(i));
            c+=child(i)->L->cols();
        }
    }
}

std::string OperatorSkeletonized::strNode() const{
    Str res("","");
    if(L!=0) res=res+"L("+((int)L->rows())+","+((int)L->cols())+") ";
    if(D!=0) res=res+"D("+((int)D->rows())+","+((int)D->cols())+") ";
    if(R!=0) res=res+"R("+((int)R->rows())+","+((int)R->cols())+") ";
    if(S!=0) res=res+"S("+((int)S->rows())+","+((int)S->cols())+") ";
    if(temporaryStorage!=0) res=res+"temp("+temporaryStorage->size()+")";
    return res;
}

/*
 * To be removed
 */
void show(MatrixXcd mat){
    std::cout<<std::endl;
    for(int i=0;i<mat.rows();i++){
        for(int j=0;j<mat.cols();j++){
            if(std::abs(mat(i,j))<1.e-10) std::cout<<".";
            else if(std::abs(mat(i,j))<1.e-2) std::cout<<"x";
            else std::cout<<"X";
        }
        
        std::cout<<std::endl;
    }
}

void show(Map<MatrixXcd> mat){
    std::cout<<std::endl;
    for(int i=0;i<mat.rows();i++){
        for(int j=0;j<mat.cols();j++){
            if(std::abs(mat(i,j))<1.e-10) std::cout<<".";
            else if(std::abs(mat(i,j))<1.e-2) std::cout<<"x";
            else std::cout<<"X";
        }
        
        std::cout<<std::endl;
    }
}


std::vector<std::complex<double> > OperatorSkeletonized::compressing;
int OperatorSkeletonized::compressingCols=0;
int OperatorSkeletonized::compressingRows=0;

OperatorSkeletonized* OperatorSkeletonized::setup(const Index* index){
    /*
     * Since currently only matching the lowest level with coefficients floor is supported
     */
    if(index->depthInFloor()!=Index::npos && index->depthInFloor()!=0) return 0;
    
    OperatorSkeletonized* result = new OperatorSkeletonized("",index);
    
    for(unsigned int i=0; i<index->childSize(); i++){
        OperatorSkeletonized* child=setup(index->child(i));
        if(child!=0) result->childAdd(child);
    }
    
    return result;
}

void OperatorSkeletonized::skeletonizeAtLevel(OperatorSkeletonized* root, int level){
    std::vector<const OperatorSkeletonized*> childrenAtLevel_;
    root->childrenAtDepth(level, childrenAtLevel_);
    std::vector<OperatorSkeletonized*> childrenAtLevel;
    for(unsigned int i=0; i<childrenAtLevel_.size(); i++) childrenAtLevel.push_back(const_cast<OperatorSkeletonized*>(childrenAtLevel_[i]));
    
    std::cout<<"SKELETONIZING LEVEL "<<level<<std::endl;
    
    std::vector<int> rowSizes, colSizes;
    for(unsigned int i=0; i<childrenAtLevel.size(); i++){
        if(childrenAtLevel[i]->childSize()==0){
            rowSizes.push_back(childrenAtLevel[i]->idx->sizeCompute());
            colSizes.push_back(rowSizes[i]);
        }else{
            int rowSize=0;
            int colSize=0;
            
            for(unsigned int j=0; j<childrenAtLevel[i]->childSize(); j++){
                rowSize+=childrenAtLevel[i]->child(j)->L->cols();
                colSize+=childrenAtLevel[i]->child(j)->R->rows();
            }
            
            rowSizes.push_back(rowSize);
            colSizes.push_back(colSize);
            
        }
        
        int tempSize=colSizes[i];
        if(rowSizes[i]>tempSize) tempSize=rowSizes[i];
        
        childrenAtLevel[i]->temporaryStorage = new std::vector<std::complex<double> >();
        childrenAtLevel[i]->temporaryStorage->resize(tempSize);
    }
    
    
    Map<MatrixXcd> matrix(compressing.data(),compressingRows,compressingCols);
    
    unsigned int col=0;
    unsigned int row=0;
    for(unsigned int i=0; i<childrenAtLevel.size(); i++){
        
        //TODO: There is a better way!
        childrenAtLevel[i]->D=new MatrixXcd(rowSizes[i],colSizes[i]);
        for(int x=0;x<rowSizes[i];x++){
            for(int y=0;y<colSizes[i];y++){
                (*childrenAtLevel[i]->D)(x,y)=matrix(row+x,col+y);
            }
        }
        
        matrix.block(row,col,rowSizes[i],colSizes[i]).setZero();
        
        row+=rowSizes[i];
        col+=colSizes[i];
    }
    
    
    std::cout<<"COMPRESSING ROW SPACE"<<std::endl;
    
    /*
     * Compress row space
     */
    Eigen::JacobiSVD<MatrixXcd> svd;
    
    MatrixXcd rowCompressed(compressingRows,compressingCols);
    compressingRows=0;
    for(unsigned int i=0; i<childrenAtLevel.size(); i++){
        /*
         * Setup matrix and calculate SVD
         */
        MatrixXcd temp(rowSizes[i],matrix.cols()-colSizes[i]);
        unsigned int c=0;
        unsigned int cTemp=0;
        for(unsigned int j=0; j<childrenAtLevel.size(); j++){
            if(i!=j){
                temp.block(0,cTemp,rowSizes[i],colSizes[j])=matrix.block(compressingRows,c,rowSizes[i],colSizes[j]);
                cTemp+=colSizes[j];
            }
            c+=colSizes[j];
        }
        
        svd.compute(temp, Eigen::ComputeThinU | Eigen::ComputeThinV);
        MatrixXcd L = svd.matrixU();
        MatrixXcd Sblock = svd.matrixV().adjoint();
        VectorXd singularValues = svd.singularValues();
        
        /*
         * Lowrank approximate
         */
        unsigned int rank=0;
        for(; rank<singularValues.size(); rank++){
            if(singularValues(rank)<1. && rank>0) break;
            
            for(unsigned int n=0; n<Sblock.cols(); n++){
                Sblock(rank,n)*=singularValues(rank);
            }
        }
        
        
        /*
         * Store L
         */
        //TODO: There is a better way
        childrenAtLevel[i]->L = new MatrixXcd(L.rows(),rank);
        for(unsigned int x=0; x<L.rows(); x++){
            for(unsigned int y=0; y<rank; y++){
                (*childrenAtLevel[i]->L)(x,y)=L(x,y);
            }
        }
        
        /*
         * Store Sblock
         */
        c=0;
        cTemp=0;
        for(unsigned int j=0; j<childrenAtLevel.size(); j++){
            if(i!=j){
                rowCompressed.block(compressingRows,c,rank,colSizes[j])=Sblock.block(0,cTemp,rank,colSizes[j]);
                cTemp+=colSizes[j];
            }
            c+=colSizes[j];
        }
        
        rowSizes[i]=rank;
        compressingRows+=rank;
    }
    
    std::cout<<"COMPRESSING COL SPACE"<<std::endl;
    
    
    /*
     * Compress col space
     */
    compressing.resize(compressingRows*compressingCols);
    Map<MatrixXcd> _compressing(compressing.data(), compressingRows, compressingCols);
    compressingCols=0;
    for(unsigned int i=0; i<childrenAtLevel.size(); i++){
        /*
         * Setup matrix and calculate SVD
         */
        MatrixXcd temp(matrix.rows()-rowSizes[i],colSizes[i]);
        unsigned int c=0;
        unsigned int cTemp=0;
        for(unsigned int j=0; j<childrenAtLevel.size(); j++){
            if(i!=j){
                temp.block(cTemp,0,rowSizes[j],colSizes[i])=rowCompressed.block(c,compressingCols,rowSizes[j],colSizes[i]);
                cTemp+=rowSizes[j];
            }
            c+=rowSizes[j];
        }
        
        svd.compute(temp, Eigen::ComputeThinU | Eigen::ComputeThinV);
        MatrixXcd Sblock = svd.matrixU();
        MatrixXcd R = svd.matrixV().adjoint();
        VectorXd singularValues = svd.singularValues();
        
        /*
         * Lowrank approximate
         */
        unsigned int rank=0;
        for(; rank<singularValues.size(); rank++){
            if(singularValues(rank)<1. && rank>0) break;
            
            for(unsigned int n=0; n<Sblock.rows(); n++){
                Sblock(n,rank)*=singularValues(rank);
            }
        }
        
        
        std::cout<<Sblock.cols()<<"->"<<rank<<std::endl;
        
        /*
         * Store R
         */
        //TODO: There is a better way
        childrenAtLevel[i]->R = new MatrixXcd(rank,R.cols());
        for(unsigned int x=0; x<rank; x++){
            for(unsigned int y=0; y<R.cols(); y++){
                (*childrenAtLevel[i]->R)(x,y)=R(x,y);
            }
        }
        
        /*
         * Store Sblock
         */
        c=0;
        cTemp=0;
        for(unsigned int j=0; j<childrenAtLevel.size(); j++){
            if(i!=j){
                _compressing.block(c,compressingCols,rowSizes[j],rank)=Sblock.block(cTemp,0,rowSizes[j],rank);
                cTemp+=rowSizes[j];
            }
            c+=rowSizes[j];
        }
        
        colSizes[i]=rank;
        compressingCols+=rank;
        
    }
    
    std::cout<<"Skeletonized"<<std::endl;
    show(_compressing);
    
    //TODO: Does this work? (column major vs row major)
    compressing.resize(compressingRows*compressingCols);
    
    std::cout<<"DONE"<<std::endl;
}


OperatorSkeletonized* OperatorSkeletonized::skeletonize(OperatorAbstract* base){
    if(base->iIndex!=base->jIndex) ABORT("Skeletonization currently only supported on equal lhs and rhs indices");
    
    OperatorSkeletonized* result = setup(base->iIndex);
    
    int dim = base->iIndex->sizeCompute();
    
    base->matrix(compressing);
    
    //MatrixXcd random = MatrixXcd::Random(dim,dim);
    
    MatrixXcd check(dim,dim);
    for(unsigned int i=0; i<dim; i++){
        for(unsigned int j=0; j<dim; j++){
            check(i,j)=compressing[i*dim+j];
            //compressing[i*dim+j]=random(i,j);
            //check(i,j)=random(i,j);
        }
    }
    
    compressingRows = dim;
    compressingCols = dim;
    
    /*
     * Skeletonize matrix starting from lowest level
     */
    OperatorSkeletonized* tmp= result;
    while(tmp->childSize()!=0)tmp=tmp->descend(1);
    unsigned int depth = tmp->depth();
    for(unsigned int i=depth;i>0;i--){
        skeletonizeAtLevel(result, i);
    }
    
    /*
     * Store skeletonized matrix in root
     */
    
    //TODO: There is a better way!
    Map<MatrixXcd> compressed(compressing.data(), compressingRows, compressingCols);
    result->S=new MatrixXcd(compressingRows, compressingCols);
    
    for(int x=0;x<compressingRows;x++){
        for(int y=0;y<compressingCols;y++){
            (*result->S)(x,y)=compressed(x,y);
        }
    }
    
    
    /*
     * Additional: Checks etc
     */
    std::cout<<result->str()<<std::endl;
    
    std::cout<<"CHECK... ";
    std::vector<std::complex<double> > checkResult;
    result->matrix(checkResult);
    
    for(unsigned int i=0; i<dim; i++){
        for(unsigned int j=0; j<dim; j++){
            if(std::abs(check(i,j)-checkResult[i*dim+j])>1.e-2) std::cout<<"ERROR "<<i<<" "<<j<<": "
                <<check(i,j)<<" "<<checkResult[i*dim+j]<<std::endl;
        }
    }
    std::cout<<"DONE ";
    
    //TODO: Free compressed
        
    return result;
}
