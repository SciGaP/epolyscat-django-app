// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORPERMUTE_H
#define OPERATORPERMUTE_H

#include <complex>
#include "operatorTree.h"

class Coefficients;
class OperatorPermute : public OperatorTree
{

    std::vector<int> _perm;
    OperatorPermute(std::string Name):OperatorTree(Name,0,0){}
public:
    /// permute around the JLevel of JIndex
    ///
    /// creates iIndex as a tree with iIndex.depth()=2, with all bases type BasisVector
    ///
    /// iIndex.childSize() is equal to number of Index's on JLevel, i.e. == JIndex.descend(JLevel).levelSize();
    /// iIndex.child(i).size() is the number of JIndex children with size at least i+1
    ///
    /// let iC(iIndex),jC(jIndex) be Coefficients and jCn
    ///
    /// the permutation is :iC.child(i).data()[j]=(j'th JIndex node at level JLevel).data()[iPos]
    /// at present, the code stops unless both, iC and jC, are stored continguously
    OperatorPermute(std::string Name, size_t JLevel, const Index* JIndex, bool MapFromIndex=true);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const override;
    std::string strNode(int Digits=-1) const override;
};

#endif // OPERATORPERMUTE_H
