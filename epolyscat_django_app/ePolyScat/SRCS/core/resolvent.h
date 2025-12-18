// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef RESOLVENT_H
#define RESOLVENT_H

#include <complex>

#include "operatorAbstract.h"
#include "tree.h"
#include "qtEigenSparse.h"

class TriFactor;
class OperatorTree;
class Resolvent : public Tree<Resolvent>, public OperatorAbstract
{
    std::complex<double> _z;
    TriFactor * factor; // the actual factorization at the leaf level
//    Eigen::
//    * SparseLU<SparseMatrix<scalar, ColMajor>, COLAMDOrdering<Index> >   solver;
//    * // fill A and b;
//    * // Compute the ordering permutation vector from the structural pattern of A
//    * solver.analyzePattern(A);
//    * // Compute the numerical factorization
//    * solver.factorize(A);
//    * //Use the factors to solve the linear system
//    * x = solver.solve(b);
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> > > solver;
public:
    ~Resolvent();
    Resolvent():factor(0){}
    /// Resolvent of Op at Z, i.e. (Op-Z)^-1
    Resolvent(const OperatorTree * Op, std::complex<double> Z);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

    std::complex<double> z() const {return _z;}

    bool verify(const OperatorTree * Op, double Epsilon=1.e-12) const;
};

#endif // RESOLVENT_H
