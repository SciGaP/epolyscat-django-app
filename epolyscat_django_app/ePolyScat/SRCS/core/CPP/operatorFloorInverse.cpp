// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorFloorInverse.h"
#include "index.h"
#ifdef _HAS_LAPACKE_
#include "lapacke.h"
#endif

OperatorFloorInverse::OperatorFloorInverse(const Index *Idx, unsigned int SubD, unsigned int SuperD, bool BandOvr) : OperatorFloor("FloorInverse"){
    _rows = Idx->sizeCompute();
    _cols = Idx->sizeCompute();
    bandOvr = BandOvr;
    UseMatrix mat;
    Idx->overlap()->matrix(mat);
    ipiv.resize(std::min(mat.rows(),mat.cols()));
#ifdef _HAS_LAPACKE_
    if(bandOvr){
        subD = SubD; superD = SuperD;
        lumat = new UseMatrix(mat.reband(subD,superD+subD,true));
        LAPACKE_zgbtrf(LAPACK_COL_MAJOR,mat.rows(),mat.cols(),lumat->subD(),lumat->superD()-lumat->subD(),
                                lumat->data(),lumat->leadDim(),ipiv.data());
    }
    else{
        lumat = new UseMatrix(mat);
        LAPACKE_zgetrf(LAPACK_COL_MAJOR,lumat->rows(),lumat->cols(),lumat->data(),lumat->leadDim(),ipiv.data());
    }
#else
    _lu.reset(new Eigen::FullPivLU<Eigen::MatrixXcd>(Idx->overlap()->matrix()));
#endif

}

OperatorFloorInverse::~OperatorFloorInverse(){
    if(lumat!=0) delete lumat;
}

void OperatorFloorInverse::axpy(const std::complex<double> & Alfa, const std::complex<double> *X, unsigned int SizX,
                  const std::complex<double> & Beta, std::complex<double> *Y, unsigned int SizY) const{
    if(Alfa==0.){scale(Beta,Y,SizY); return;}
    else if(Alfa!=1.) ABORT("Not yet implemented for Alfa!=1, except Alfa==0");
    if(SizX!=SizY) ABORT("SizX != SizY. Cannot use axpy.");

    std::vector<std::complex<double> > Z;
    for(unsigned int i=0;i<SizX;i++) Z.push_back(X[i]);

#ifdef _HAS_LAPACKE_
    const char trans = 'N';
    if(bandOvr){
        LAPACKE_zgbtrs(LAPACK_COL_MAJOR,trans,lumat->cols(),lumat->subD(),lumat->superD()-lumat->subD(),1,
                                lumat->data(),lumat->leadDim(),ipiv.data(),
                                Z.data(),SizX);
    }
    else{
        LAPACKE_zgetrs(LAPACK_COL_MAJOR,trans,lumat->cols(),1,lumat->data(),lumat->leadDim(),ipiv.data(),
                                Z.data(),SizX);
    }
#else
    _lu->solve(Eigen::Map<Eigen::MatrixXcd>(Z.data(),SizY,1));
#endif

    if(Beta==0.) for(unsigned int i=0;i<SizY;i++) Y[i] = Z[i];
    else for(unsigned int i=0;i<SizY;i++) Y[i] = Z[i] + Beta * Y[i];
}
