// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
//
//  operatorSkeletonized.h
//  tRecX
//
//  Created by Jonas Bucher on 19.04.17.
//  Copyright © 2017 youdontneedtoknow. All rights reserved.
//

#ifndef OPERATOR_SKELETONIZED_H
#define OPERATOR_SKELETONIZED_H

#include <stdio.h>
#include "operatorAbstract.h"
#include "tree.h"
#include "qtEigenDense.h"

/**
 * Warning! Abandoned
 */
class OperatorSkeletonized: public OperatorAbstract, public Tree<OperatorSkeletonized>{
    const Index* idx;
    
    /*
     * Either these three matrices are present marking the node as not the root
     */
    Eigen::MatrixXcd* R;
    Eigen::MatrixXcd* L;
    Eigen::MatrixXcd* D;
    
    //STYLE remove if used only for construction
    std::vector<std::complex<double> >* temporaryStorage;
    
    /*
     * Or this matrix is present marking the node as the root
     */
    Eigen::MatrixXcd* S;
    
    static std::vector<std::complex<double> > compressing;
    static int compressingCols;
    static int compressingRows;
    
    
public:
    OperatorSkeletonized(): R(0), L(0), D(0), temporaryStorage(0), S(0) {}
    OperatorSkeletonized(std::string name, const Index* index): OperatorAbstract(name, index, index), idx(index), L(0), R(0), D(0), S(0), temporaryStorage(0){}
    
    void apply(std::complex<double> alpha, const Coefficients &x, std::complex<double> beta, Coefficients &y) const;
    void compressStoringDiagonals(std::complex<double>* target, const Coefficients& x) const;
    void decompressRestoringDiagonals(std::complex<double>* source, Coefficients& y) const;
    
    std::string str() const {return Tree<OperatorSkeletonized>::str();}
    std::string strNode() const;
    
    static OperatorSkeletonized* skeletonize(OperatorAbstract* base);
    
private:
    static OperatorSkeletonized* setup(const Index* index);
    static void skeletonizeAtLevel(OperatorSkeletonized* root, int level);
    
};

#endif
