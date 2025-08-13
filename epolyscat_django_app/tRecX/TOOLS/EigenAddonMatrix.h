// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef EIGENADDONMATRIX_H
#define EIGENADDONMATRIX_H

/// WARNING: conflicts some headers of stdlib, eg. with <vector> in gcc-7, <stdio> in gcc 4.8.5
///
/// Problems frequently appear with "namespace std _GLIBCXX_VISIBILITY(default)"
/// !!! DO NOT INCLUDE THESE HEADERS HERE !!!

/// Usage: include (in the following sequence!)
/// #define EIGEN_MATRIX_PLUGIN "EigenAddonMatrix.h"
/// #include "qtEigenDense.h"

/// set near-zeros =0, near-zero: < max(Eps*rowNorm,Eps*colNorm,EpsAbs)
inline Matrix<Scalar,-1,-1> & purge(double Eps=1.e-12,double EpsAbs=1.e-14){
    if(Eps<=0. and EpsAbs<=0.)return *this; // zero threshold, no purge

    //HACK: as of gcc-6, plugins create ugly conflicts with std::vector - use array instead
    // std::vector<double> epsI(this->rows(),0.),epsJ(this->cols(),0.);
    double epsI[this->rows()],epsJ[this->cols()];
    for(int k=0;k<this->rows();k++)epsI[k]=std::max(this->row(k).norm()*Eps,EpsAbs);
    for(int k=0;k<this->cols();k++)epsJ[k]=std::max(this->col(k).norm()*Eps,EpsAbs);
    for(unsigned int j=0;j<this->cols();j++) {
        for(unsigned int i=0;i<this->rows();i++){
            if(std::abs(this->operator()(i,j).imag())<std::max(epsI[i],epsJ[j]))
                this->operator()(i,j)=std::complex<double>(this->operator()(i,j).real(),0.);
            if(std::abs(this->operator()(i,j).real())<std::max(epsI[i],epsJ[j]))
                this->operator()(i,j)=std::complex<double>(0.,this->operator()(i,j).imag());
        }
    }
    return *this;
}

#endif // EIGENADDONMATRIX_H
