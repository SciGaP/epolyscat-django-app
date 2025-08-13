// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORDIAGSPERMUTED_H
#define OPERATORDIAGSPERMUTED_H

#include "operatorTree.h"

/// Operator compose of permutations and a block-diagonal operator
///
/// Y = Op X = Perm^-1 OpBlockDiag Perm X
class OperatorDiagsPermuted: public OperatorTree
{
    //HACK
    friend class OperatorRALL;

    std::unique_ptr<Coefficients> _rhs,_lhs;
    OperatorTree *perm,*back,*diag;
    void construct(std::map<const Index*, std::map<const Index*, Eigen::VectorXcd> > &Diags);

    // suppress copy and assignment
    OperatorDiagsPermuted(const OperatorDiagsPermuted&){};
    OperatorDiagsPermuted& operator=(const OperatorDiagsPermuted&){return *this;};

    /// (DEBUG only) from Operator tree that is composed of diagonal blocks
    OperatorDiagsPermuted(OperatorTree& Op, double DiagonalThreshold /** maximal off-diagonal elements */);
public:
    ~OperatorDiagsPermuted();
    /// construct from Index-labeled matrix of diagonal blocks
    OperatorDiagsPermuted(std::string Name, const Index *IIndex, const Index*JIndex,
                       std::map<const Index*,std::map<const Index*,Eigen::VectorXcd>> &Diags);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    bool isEfficient() const{return true;}
    bool symmetrize(); ///< checks for symmetry and symmetrizes explicitly
};


#endif // OPERATORDIAGSPERMUTED_H
