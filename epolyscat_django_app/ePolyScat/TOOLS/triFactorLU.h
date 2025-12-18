// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef LU_H
#define LU_H
#include <memory>
#include "triFactor.h"
#include "abort.h"
#include <Sparse>

/// \ingroup Linalg
/// \brief Abstract for LU matrix factorization
class TriFactorLU: public TriFactor
{
    TriFactorLU(const TriFactorLU & other){ABORT("do not use copy constructor");}
    TriFactorLU(const UseMatrix &M, bool ImproveIteratively=false);
    std::shared_ptr<Eigen::FullPivLU<Eigen::MatrixXcd>> _lu;
    std::shared_ptr<Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>> _slu;
public:
    ~TriFactorLU();
    TriFactorLU(){shapeM=0;}
    UseMatrix & solve(const char Trans, UseMatrix &Rhs) const;
    UseMatrix inverse() const;
    std::complex<double> det() const;

    void reFactor(const UseMatrix & M);
};

#endif // LU_H
