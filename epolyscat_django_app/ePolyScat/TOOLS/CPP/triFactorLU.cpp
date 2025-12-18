// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "triFactorLU.h"
#include "useMatrix.h"
#include "abort.h"
#ifdef _HAS_LAPACKE_
#include "lapacke.h"
#endif
#include "tools.h"
#include "str.h"
#include <Sparse>

using namespace std;

#ifdef _HAS_LAPACKE_
static std::string lapackInfo(std::string Routine, int info){
    return Routine+": info= "+tools::str(info);
}
#endif

TriFactorLU::~TriFactorLU(){
    delete shapeM;
    shapeM=0;
}

TriFactorLU::TriFactorLU(const UseMatrix &M, bool ImproveIteratively)
{
    if(M.rdata!=0)ABORT("not implmenented for real yet");
#ifndef _HAS_LAPACKE_
    Eigen::SparseMatrix<std::complex<double>> sM;
    sM=Eigen::Map<Eigen::MatrixXcd>(M.data(),M.rows(),M.cols()).sparseView();
    _slu.reset(new Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>());
    _slu->analyzePattern(sM);
    _slu->factorize(sM);

//    _lu.reset(new Eigen::FullPivLU<Eigen::MatrixXcd>(Eigen::Map<Eigen::MatrixXcd>(M.data(),M.rows(),M.cols())));
#else
    shapeM=new UseMatrix::Shape(M.shape);
    pivot.resize(min(M.rows(),M.cols()));

    if(M.size()==0){
        ABORT("cannot LU decompose 0-size matrix");
    }
    else if(M.size()==1){
        // Lapack cannot handle size=1 case
        if(M(0,0).complex()==0.)ABORT("1x1 zero matrix");
        storeLU=new UseMatrix(M);
        storeLU->operator()(0,0)=1./M(0,0).complex();
        return;
    }

    unsigned int nonZero,trueSub,trueSuper;
    char dataType;
    M.diagnose(0,nonZero,trueSub,trueSuper,dataType);
    if(M.rows()<10 or (trueSub+trueSuper)>max(M.rows(),M.cols())){
        if(ImproveIteratively)mat=new UseMatrix(M);
        if(shapeM->locationIJ!=UseMatrix::full_normal)
            ABORT("DEVELOPER: full matrix - need to rearrange from banded storage");

        storeLU=new UseMatrix(M);
        int info=LAPACKE_zgetrf(LAPACK_COL_MAJOR,
                                shapeM->nrows,shapeM->ncols,storeLU->data(),shapeM->nrows,pivot.data());
        if(info!=0)ABORT(lapackInfo("zgetrf",info));
    }
    else {
        storeLU=new UseMatrix(M.reband(trueSub,trueSuper+trueSub,true));

        // I was being smart: reband does not always return a banded matrix
        int info;
        if(storeLU->isFull()){
            if(ImproveIteratively)mat=new UseMatrix(M);
            info=LAPACKE_zgetrf(LAPACK_COL_MAJOR,
                                shapeM->nrows,shapeM->ncols,storeLU->data(),shapeM->nrows,pivot.data());
        }
        else{
            if(ImproveIteratively)mat=new UseMatrix(M.reband(trueSub,trueSuper,true));
            info=LAPACKE_zgbtrf(LAPACK_COL_MAJOR,
                                M.rows(),M.cols(),trueSub,trueSuper,
                                storeLU->data(),storeLU->leadDim(),pivot.data());
        }
        if(info!=0){
            ABORT(lapackInfo("zgbtrf",info));
        }
    }
#endif
}

UseMatrix & TriFactorLU::solve(const char Trans, UseMatrix &Rhs) const{
    if(shapeM->nrows!=shapeM->ncols)ABORT("non-square matrix cannot be inverted");

#ifndef _HAS_LAPACKE_
    _slu->solve(Eigen::Map<Eigen::MatrixXcd>(Rhs.data(),Rhs.rows(),Rhs.cols()));
    return Rhs;
#else

    // save Rhs for case of iterative improvement
    UseMatrix SaveRhs;
    if(mat!=0)SaveRhs=Rhs;

    if(storeLU->size()==1){
        Rhs(0,0)*=storeLU->operator()(0,0).complex();
    }
    else if(storeLU->isFull()) {
        // full matrix
        int info=LAPACKE_zgetrs(LAPACK_COL_MAJOR,
                                Trans,shapeM->ncols,Rhs.cols(),storeLU->data(),shapeM->nrows,pivot.data(),
                                Rhs.cdata,Rhs.leadDim());
        if(info!=0)ABORT(lapackInfo("zgetrs",info));
    }
    else if (storeLU->isBand()) {
        // banded matrix
        //        lapack_int LAPACKE_zgbtrs( int matrix_layout, char trans, lapack_int n,
        //                                   lapack_int kl, lapack_int ku, lapack_int nrhs,
        //                                   const lapack_complex_double* ab, lapack_int ldab,
        //                                   const lapack_int* ipiv, lapack_complex_double* b,
        //                                   lapack_int ldb )
        int info=LAPACKE_zgbtrs(LAPACK_COL_MAJOR,
                                Trans,storeLU->cols(),storeLU->subD(),storeLU->superD()-storeLU->subD(),Rhs.cols(),
                                storeLU->data(),storeLU->leadDim(),pivot.data(),
                                Rhs.cdata,Rhs.leadDim());
        if(mat!=0){
            UseMatrix Old(Rhs);
            vector<double> ferr(Rhs.cols()),berr(Rhs.cols());
            info=LAPACKE_zgbrfs(LAPACK_COL_MAJOR,
                                'n',mat->cols(),mat->subD(),mat->superD(),Rhs.cols(),
                                mat->data(),mat->leadDim(),storeLU->data(),storeLU->leadDim(),
                                pivot.data(),SaveRhs.data(),SaveRhs.leadDim(),
                                Rhs.data(),Rhs.leadDim(),ferr.data(),berr.data());

            Str("improving: errors and change")+ferr+berr+(Old-Rhs).maxAbsVal()+Str::print;
            if(info!=0)ABORT(lapackInfo("zgbrfs",info));
        }

    }
    else
        ABORT("storage format not covered");
    return Rhs;
#endif
}

UseMatrix TriFactorLU::inverse() const{
#ifndef _HAS_LAPACKE_
    UseMatrix res=UseMatrix::Identity(_slu->rows(),_slu->cols());
    Eigen::Map<Eigen::MatrixXcd>(res.data(),res.rows(),res.cols())
            =_slu->solve(Eigen::Map<Eigen::MatrixXcd>(res.data(),res.rows(),res.cols()));
    return res;
#else
    if(shapeM->nrows!=shapeM->ncols)ABORT("non-square matrix does not have an inverse");
    UseMatrix inv;
    if(storeLU->isFull()){
        inv=UseMatrix(shapeM->nrows,shapeM->ncols);
        complex<double> *cd=inv.cdata;
        for(unsigned int i=0;i<storeLU->size();i++,cd++)*cd=storeLU->data()[i];
        int info=LAPACKE_zgetri(LAPACK_COL_MAJOR,
                                inv.cols(),inv.cdata,inv.leadDim(),pivot.data());
        if(info!=0)ABORT(lapackInfo("zgetri",info));
    }
    else {
        // there does not seem to exist LAPACK zgbtri (probably does not make much sense).
        inv=UseMatrix::Identity(storeLU->rows(),storeLU->cols());
        solve('n',inv);
    }
    return inv;
#endif
}
complex<double> TriFactorLU::det() const {
    if(shapeM->nrows!=shapeM->ncols)ABORT("non-square matrix does not have a determiant");
    ABORT("determinant not implemented for shape="+shapeM->strLocIJ());
    return complex<double>(1.);
}

void TriFactorLU::reFactor(const UseMatrix &M){
    TriFactorLU lu(M);
    swap(storeLU,lu.storeLU);
    swap(pivot,lu.pivot);
    swap(_lu,lu._lu);
    swap(_slu,lu._slu);
    delete shapeM;
    shapeM = new UseMatrix::Shape(M.shape);
}
