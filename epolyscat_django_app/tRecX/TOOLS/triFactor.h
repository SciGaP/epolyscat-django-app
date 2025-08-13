// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef TRIFACTOR_H
#define TRIFACTOR_H
#include "useMatrix.h"

/// \ingroup Linalg
/// \brief base class for triangular factorizations
class TriFactor
{
    friend class UseMatrix;
public:
    virtual ~TriFactor(){delete storeLU;delete mat;}
    TriFactor():shapeM(0),storeLU(0),mat(0){}
    virtual UseMatrix inverse() const =0;
    /// overwrite Rhs with solution, also return reference to Rhs
    virtual UseMatrix & solve(const char Trans, UseMatrix & Rhs) const =0;

    /// LU factorization allows iterative improvement of given solution
    virtual void improve(const UseMatrix & M, const UseMatrix & Rhs, UseMatrix & Solution) const{}

    virtual std::complex<double> det() const =0;
    virtual void reFactor(const UseMatrix & M)=0; /// (re-)compute the factorization of M
    void reset(); /// reset (when main matrix is modifield
    bool verify(const UseMatrix & M);
protected:
    std::vector<int> pivot;
    UseMatrix::Shape * shapeM; /// pointer to the shape of the original matrix
    UseMatrix * storeLU;
    UseMatrix * mat; /// original matrix
};

#endif // TRIFACTOR_H
