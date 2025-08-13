// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
//
// Created by Jonas Bucher on 05.04.17.
//

#ifndef OPERATORSVD_H
#define OPERATORSVD_H

#include <complex>
#include <string>

#include "operatorFloor.h"
#include "operatorTree.h"
#include "index.h"

#define _OPERATOR_SVD_OUTPUT_

/** \ingroup Structures */
/// @brief Lowrank optimized OperatorAbstract
///
/// Decomposes a given OperatorAbstract into U * diag * V^\dagger, U, V unitary, diag diagonal.
/// From this an effective rank is determined based on the singular values in diag giving an optimized
/// matrix-vector operations count
class OperatorSVD: public OperatorFloor{
private:
    const static double SINGULAR_VALUE_CUTOFF; ///< Singular values below this level will be removed
    
    const Index* iIndex;
    const Index* jIndex;

    std::vector<std::complex<double> >* U;
    std::vector<std::complex<double> >* V;
    std::vector<std::complex<double> >* diag;

    // Information about the approximation, can be extended later to not just use a constant number of svs
    unsigned int rank;

    void check(OperatorAbstract* base);

protected:
    void axpy(const std::complex<double> & Alfa, const std::complex<double>*X, unsigned int SizX,
              const std::complex<double> & Beta,       std::complex<double>*Y, unsigned int SizY) const;
    void axpy(std::complex<double> Alfa, const std::vector<std::complex<double> >& X,
              std::complex<double> Beta,       std::vector<std::complex<double> >& Y) const;
public:
    OperatorSVD(OperatorAbstract* base, std::string Definition=""); ///< SVD decomposes base
    virtual ~OperatorSVD(){ delete U; delete V; delete diag; }

    std::vector<double> getSingularValues() const; ///< Getter

    // TODO
    void pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const;//{ ABORT("Not implemented"); }

    long applyCount() const{ return (_rows+_cols)*rank + diag->size(); }
     
    std::string strInfo() const;
    
    /// @brief Optimize OperatorTree at each level
    ///
    /// The operator tree is mutated by this method. Starting from the top, SVD decompositions
    /// of the optree nodes are compared in performance to their original counterparts and placed into the tree if faster. All determined singuar 
    /// values are written into output (if output is given).
    static void optimize(OperatorTree* base, 
                         bool respectFloorLevel=false, ///< If true, only operator floors will be optimized 
                         std::ostream* output=nullptr);
    
    /// @brief Simple test of the class
    static void test();

    /// @brief Write singular values to ostream
    void write(std::ostream* output) const;
    
    
    /*
     ******************************* OLD CODE ************************************
     */
    void apply(std::complex<double> alpha, const Coefficients &x, std::complex<double> beta, Coefficients &y) const;
    std::string strNode() const;
};

#endif //TRECX_OPERATORSVD_H
