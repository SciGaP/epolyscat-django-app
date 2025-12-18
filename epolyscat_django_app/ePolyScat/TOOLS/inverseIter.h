// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INVERSEITER_H
#define INVERSEITER_H
#include <vector>
#include "useMatrix.h"
#include "arpack.h"

class InverseIter:public Arpack
{
    UseMatrix aMinusEgM;
    const UseMatrix &mat,&ovr;
    UseMatrix yTemp,xTemp;

    void apply(const std::complex<double> *X, std::complex<double> *Y);

public:
    static void test();

    InverseIter(const UseMatrix &Mat, const UseMatrix & Ovr)
        :Arpack(1.e-12,100),mat(Mat),ovr(Ovr),xTemp(UseMatrix(Mat.cols(),1)),yTemp(UseMatrix(Mat.cols(),1)){_lvec=yTemp.size();}

    /// \brief eigen - find several eigenvalues:  Mat Evec[k] = Ovr Evec[v] Eval[k]
    /// \param Eval  - initial guess valeus, Eval.size()  > 0
    /// \param Evec  - columns are eigenvectors, non-empty: guess vectors
    void eigen(UseMatrix & Eval, UseMatrix & Evec);

};

#endif // INVERSEITER_H
