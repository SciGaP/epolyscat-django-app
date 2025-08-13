// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DERIVATIVEBLOCK_H
#define DERIVATIVEBLOCK_H

#include <vector>
#include <complex>


class Coefficients;
class OperatorFloor;
class ParallelFloor;
class OperatorTree;
class DerivativeBlock {
public:
    static bool lessEqual(const DerivativeBlock & A, const DerivativeBlock & B);
    static bool less(const DerivativeBlock & A, const DerivativeBlock & B);
    DerivativeBlock(std::vector<unsigned int> BlockSort,
                    const OperatorTree* OLeaf, Coefficients * ICoef, Coefficients*JCoef, double ONorm, double *XNorm);

    std::vector<int> _info; // floor structural info
    const OperatorTree * oLeaf;
    std::vector<Coefficients*> cInOut;
#ifdef _OPENMP
    mutable std::complex<double>* _alt; // temporary storag in OMP loop
    mutable bool _skip;                 // skip for norm reasons
#endif
    //NOTE: factor duplicates oLeaf->floor()->factor - should be removed
    double eps;    // rejection threshold
    double* xNorm; // pointer to norm of pX

    std::vector<unsigned int> blockSort; // multi-index that will be used for sorting
    double load() const;
private:
    DerivativeBlock();
};


#endif // DERIVATIVEBLOCK_H
