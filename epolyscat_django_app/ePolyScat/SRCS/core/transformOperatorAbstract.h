// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef TRANSFORM_OPERATOR_ABSTRACT_H
#define TRANSFORM_OPERATOR_ABSTRACT_H

#include <map>

#include "useMatrix.h"
#include "coefficients.h"

class OperatorAbstract;
class OperatorTree;
class Index;

class TransformOperatorAbstract{
private:
    const OperatorAbstract* op;
    const OperatorAbstract* leftTrafo;
    const OperatorAbstract* rightTrafo;
    OperatorTree* transformed;

    std::vector<int> optimize__diagonalLevels;

    void build(OperatorTree* optree, OperatorAbstract* op, UseMatrix* mat=0, const Index* matIIndex=0, const Index* matJIndex=0);
public:
    /*
     * I <--leftTrafo-- i <-----transformed------ j <--rightTrafo-- J
     * I <--------------------------op----------------------------- J
     */
    TransformOperatorAbstract(const OperatorAbstract* Op, const OperatorAbstract* LeftTrafo, const OperatorAbstract* RightTrafo);

    void optimize__addDiagonalLevel(int l){ optimize__diagonalLevels.push_back(l); }

    void transform();
    bool makeDiagonal(std::map<std::string, double>& report, double eps=1.e-12); // \tilde C= b^\dagger A^{\sim 1}b - C
    void check();
    OperatorTree* getTransformed(){ return transformed; }
    OperatorAbstract* getTransformedWithTransformations();
};







#endif //TRANSFORM_OPERATOR_ABSTRACT_H
