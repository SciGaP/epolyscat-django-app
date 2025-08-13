// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

//
// Created by Jonas Bucher on 04.04.17.
//


#include "operatorSVD.h"

#include <vector>
#include <string>
#include "tools.h"
#include <ctime>

#include "mpiWrapper.h"
#include "readInput.h"
#include "printOutput.h"
#include "timer.h"

#include "coefficients.h"
#include "index.h"
#include "basisAbstract.h"

#include "qtEigenDense.h"
#include <Core>
#ifdef _HAS_LAPACKE_
#include "lapacke.h"
#endif

//#define _OPERATOR_SVD_NEWCODE_


const double OperatorSVD::SINGULAR_VALUE_CUTOFF = 1.e-4;

TIMER(svDecompose,)
#ifdef _OPERATOR_SVD_USE_EIGEN_
static void svDecompose(std::complex<double>* A, unsigned int rows, unsigned int cols, double* s, std::complex<double>* U, std::complex<double>* V){
    STARTDEBUG(svDecompose)
            Eigen::JacobiSVD<Eigen::MatrixXcd> svd;
    
    svd.compute(Eigen::Map<Eigen::MatrixXcd>(A, rows, cols), Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    Eigen::Map<Eigen::MatrixXcd>(U, rows, rows)          = svd.matrixU();
    Eigen::Map<Eigen::MatrixXcd>(V, cols, cols)          = svd.matrixV().adjoint();
    Eigen::Map<Eigen::VectorXd> (s, std::min(rows,cols)) = svd.singularValues();
    STOPDEBUG(svDecompose)
}
#else
static void svDecompose(std::complex<double>* A, unsigned int rows, unsigned int cols, double* s, std::complex<double>* U, std::complex<double>* V){
#ifdef _HAS_LAPACKE_
    STARTDEBUG(svDecompose)
            LAPACKE_zgesdd(/* matrix order */ LAPACK_COL_MAJOR, /* jobz */ 'A', /* m */ rows,/* n */ cols, /* A */ A, /* LDA */ rows,
                           /* s */ s,/* U */ U, /* LDU */ rows, /* VT */ V, /* LDVT */ cols);
    STOPDEBUG(svDecompose);
#else
    DEVABORT("reimplement w/o LAPACKE")
#endif
}
#endif

TIMER(operatorSVD_matrix,)
void OperatorSVD::pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const{
    DEVABORT("Not implemented");
}

OperatorSVD::OperatorSVD(OperatorAbstract *base, std::string Definition): OperatorFloor(0,0,"SVD"), iIndex(base->iIndex), jIndex(base->jIndex) {  

    _rows = iIndex->sizeStored();
    _cols = jIndex->sizeStored();


    std::vector<std::complex<double> > mat;
    
    STARTDEBUG(operatorSVD_matrix)
            base->matrix(mat);
    STOPDEBUG(operatorSVD_matrix)

            this->diag=new std::vector<std::complex<double> >();
#ifndef _OPERATOR_SVD_NEWCODE_
    /*
     * Subtract diagonals
     */
    for(unsigned int i=0; i<_rows && i<_cols; i++){
        this->diag->push_back(mat[i+i*_rows]);
        mat[i+i*_rows]=0;
    }
#else
    // TODO There is a better way, since U and Vdagger are not needed at this step
    std::vector<std::complex<double> > mat1(mat);
    std::vector<double> s1(std::min(_rows, _cols));
    std::vector<std::complex<double> > U1(_rows*_rows);
    std::vector<std::complex<double> > Vdagger1(_cols*_cols);

    svDecompose(mat1.data(), _rows, _cols, s1.data(), U1.data(), Vdagger1.data());

    double lambda=0.;
    
    /*
     * Find optimal value for lambda
     */
    unsigned int cMax=0;
    unsigned int cStart=0;
    for(unsigned int i=1; i<s1.size(); i++){
        if(s1[cStart]-s1[i]<2*SINGULAR_VALUE_CUTOFF){
            if(i-cStart > cMax){
                lambda = s1[cStart]-s1[i];
                cMax = i-cStart;
            }
        }else{
            cStart=i;
        }
    }

    /*
     * Subtract lambda times identity
     */
    for(unsigned int i=0; i<_rows && i<_cols; i++){
        this->diag->push_back(lambda);
        mat[i+i*_rows]-=lambda;
    }
#endif // _OPERATOR_SVD_NEWCODE_
    
    /*
     * Decompose
     */
    std::vector<double> s(std::min(_rows, _cols));
    std::vector<std::complex<double> > U(_rows*_rows);
    std::vector<std::complex<double> > Vdagger(_cols*_cols);
    
    svDecompose(mat.data(), _rows, _cols, s.data(), U.data(), Vdagger.data());

    /*
     * Absorb singular values into U
     */
    for(int i=0;i<s.size();i++){
        for(unsigned int j=0; j<_rows; j++) U[j+i*_rows] *= s[i];
    }

    /*
     * Approximate by setting rank
     */
    rank = 0;
    for(rank=1;rank<s.size();rank++){
        if(s[rank]<SINGULAR_VALUE_CUTOFF) break;
    }
    
    /*
     * Store data
     */
    this->U=new std::vector<std::complex<double> >(_rows*_rows);
    this->V=new std::vector<std::complex<double> >(_cols*_cols);
    
    Eigen::Map<Eigen::MatrixXcd>(this->U->data(), _rows, _rows) = Eigen::Map<Eigen::MatrixXcd>(U.data(), _rows, _rows);
    Eigen::Map<Eigen::MatrixXcd>(this->V->data(), _cols, _cols) = Eigen::Map<Eigen::MatrixXcd>(Vdagger.data(), _cols, _cols).adjoint();
    
    PrintOutput::message(Str("SVD of (")+_rows+"x"+_cols+"):"+std::min(_rows, _cols)+"->"+rank);
    
    // OPTIONAL
    check(base);
}
void OperatorSVD::axpy(const std::complex<double> & Alfa, const std::complex<double>*X, unsigned int SizX,
                       const std::complex<double> & Beta,       std::complex<double>*Y, unsigned int SizY) const{
    DEVABORT("Not implemented");
}

void OperatorSVD::axpy(std::complex<double> alpha, const std::vector<std::complex<double> >& x,
                       std::complex<double> beta,        std::vector<std::complex<double> >& y) const{


    if(beta!=1.){
        scale(beta,y);
    }

    if(alpha==0.) return;

    for(unsigned int i=0; i<diag->size(); i++) y[i]+= alpha* (*diag)[i]*x[i];
    Eigen::Map<Eigen::VectorXcd>(y.data(),_rows) += alpha* (Eigen::Map<Eigen::MatrixXcd>(U->data(),_rows,rank) *
                                                            (Eigen::Map<Eigen::MatrixXcd>(V->data(),_cols,rank).adjoint() * Eigen::Map<Eigen::VectorXcd>(const_cast<std::vector<std::complex<double> >& >(x).data(), _cols)));
}


std::string OperatorSVD::strInfo() const{
    Str s("","");
    s=s+_rows+"x(SVD: "+rank+")x"+_cols;
    return s;
}


std::vector<double> OperatorSVD::getSingularValues() const {
    std::vector<double> res;

    for(unsigned int i=0; i<std::min(iIndex->sizeStored(), jIndex->sizeStored()); i++){
        double norm = 0.;
        for(unsigned int j=0; j<iIndex->sizeStored(); j++) norm+= std::pow(std::abs((*U)[j+i*iIndex->sizeStored()]),2);
        res.push_back(std::sqrt(norm));
    }

    return res;
}


TIMERRECURSIVE(optimize,)
void OperatorSVD::optimize(OperatorTree* base, bool respectFloorLevel, std::ostream* output){
    STARTDEBUG(optimize)

            bool testAtThisLevel = true;
    bool continueTraversal = true;
    
    unsigned int baseApplyCount = base->applyCount();
    
    unsigned int rows = base->iIndex->sizeStored();
    unsigned int cols = base->jIndex->sizeStored();


    if(respectFloorLevel and not (base->iIndex->hasFloor() and base->jIndex->hasFloor())){
        testAtThisLevel = false;
    }

    /*
     * Prevent SVD of too large matrices
     */
    if(rows>1000 or cols>1000){
        testAtThisLevel = false;
    }
    
    /*
     * Prevent SVD of diagonal matrices
     */
    if(baseApplyCount<=rows or baseApplyCount<=cols){
        testAtThisLevel = false;
        continueTraversal = false;
    }
    

    if(testAtThisLevel){
        OperatorSVD* opsvd = new OperatorSVD(base, base->def());
        if(output!=0) opsvd->write(output);
        
        if(opsvd->applicationCost()<base->applicationCost()){
            base->oFloor = opsvd;
            for(unsigned int i=base->childSize(); i>0; i--){
                base->childErase(i-1);
            }

            STOPDEBUG(optimize)
                    return;
        }else{
            delete opsvd;
        }
    }
    
    if(continueTraversal){
        for(unsigned int i=0; i<base->childSize(); i++){
            optimize(base->child(i), respectFloorLevel, output);
        }
    }
    
    STOPDEBUG(optimize)
}

void OperatorSVD::write(std::ostream* output) const{
    
    std::vector<const Index*> iPath=iIndex->path();
    std::vector<const Index*> jPath=jIndex->path();
    iPath.push_back(iIndex);
    jPath.push_back(jIndex);
    
    for(unsigned int i=0; i<iPath.size()-1; i++){
        if(i!=0) *output<<"-";
        
        std::string dvr="";
        DEVABORT("BasisSet cannot be used any more");
//        const BasisSet* iBas = iPath[i]->basisSet();
//        if(iBas->name().find("DVR")!=std::string::npos){
//            UseMatrix iGrid, iWeights;
//            iBas->dvrRule(iGrid, iWeights);
//            dvr=" DVR: "+tools::str(iGrid(iPath[i+1]->nSibling()).real());
//        }
        *output<<iPath[i]->axisName()<<"("<<iPath[i+1]->nSibling()<<dvr<<")";
    }
    
    *output<<",";
    
    for(unsigned int i=0; i<jPath.size()-1; i++){
        if(i!=0) *output<<"-";
        
        std::string dvr="";
        DEVABORT("BasisSet cannot be used any more");
//        const BasisSet* jBas = jPath[i]->basisSet();
//        if(jBas->name().find("DVR")!=std::string::npos){
//            UseMatrix jGrid, jWeights;
//            jBas->dvrRule(jGrid, jWeights);
//            dvr=" DVR: "+tools::str(jGrid(jPath[i+1]->nSibling()).real());
//        }
        *output<<jPath[i]->axisName()<<"("<<jPath[i+1]->nSibling()<<dvr<<")";
    }
    
    *output<<",";

    std::vector<double> singularValues = getSingularValues();
    
    for(unsigned int i=0; i<singularValues.size(); i++){
        if(i!=0) *output<<",";
        *output<<std::setprecision(10)<<singularValues[i];
    }
    
    *output<<std::endl;
}

void OperatorSVD::check(OperatorAbstract* base){
    
    Coefficients x(jIndex);
    Coefficients y(iIndex);
    Coefficients yBase(iIndex);

    y.setToZero();
    yBase.setToZero();

    x.treeOrderStorage();
    y.treeOrderStorage();
    yBase.treeOrderStorage();

    std::vector<std::complex<double> *> pX, pY, pYBase;
    x.pointerToC(pX);
    y.pointerToC(pY);
    yBase.pointerToC(pYBase);


    double maxAbsError=-1.;
    double maxRelError=-1.;

    for(unsigned int i=0; i<25; i++){
        x.setToRandom();

        apply(1.,x,1.,y);
        base->apply(1.,x,1.,yBase);

        for(unsigned int i=0; i<pY.size(); i++){
            double absError=std::abs(*pY[i]-*pYBase[i]);
            double relError=absError/std::max(std::abs(*pY[i]),std::abs(*pYBase[i]));

            if(absError>maxAbsError) maxAbsError=absError;
            if(relError>maxRelError) maxRelError=relError;
        }

    }

    PrintOutput::message(Str("SVD Test results for 25 vectors: Max errors (abs/rel)e-5: ")+1.e5*maxAbsError+"/"+1.e5*maxRelError);
}




void OperatorSVD::test(){
    // TODO!
}


/*
 ************************** OLD CODE ****************************
 */
void OperatorSVD::apply(std::complex<double> alpha, const Coefficients &x, std::complex<double> beta, Coefficients &y) const{
    if(x.idx()!=jIndex) ABORT("rhs indices don't match");
    if(y.idx()!=iIndex) ABORT("lhs indices don't match");
    
    if(beta==0.)y.setToZero();
    else if(beta!=1.)y*=beta;
    beta=1.;
    
#ifdef _DEBUG_
    std::vector<std::complex<double>* > pX,pY;
    const_cast<Coefficients&>(x).pointerToC(pX);
    y.pointerToC(pY);
    
    for(int i=0; i<pX.size()-1;i++){
        if(pX[i+1]-pX[i]!=1) ABORT("SVD only supported for contiguous storage");
    }

    for(int i=0; i<pY.size()-1;i++){
        if(pY[i+1]-pY[i]!=1) ABORT("SVD only supported for contiguous storage");
    }
#endif

    OperatorFloor::axpy(alpha,const_cast<Coefficients&>(x).data(),x.size(),beta,y.data(),y.size());
}

std::string OperatorSVD::strNode() const{
    Str s("","");
    
    s=s+"<"+(iIndex->index())+"|"+(jIndex->index())+"> ("+iIndex->sizeStored()+"x"+jIndex->sizeStored()+") SVD("+rank+")"+(iIndex->axisName());
    if(iIndex->axisName()!=jIndex->axisName())s=s+"-"+(jIndex->axisName());
    
    return s;
}











