// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATOR_TUCKER_H
#define OPERATOR_TUCKER_H

#include <vector>
#include <complex>

#include "useMatrix.h"
#include "operatorAbstract.h"
#include "operatorTree.h"

#include "operatorAbstractProduct.h"


/**
 * Warning! Abandoned
 *
 */
class OperatorTucker: public OperatorAbstractProduct{
private:
    static void ensureTensorProduct(const Index* idx, const Index* base, std::vector<int>& dims);
    static void calculateDensityMatrix(UseMatrix& mat, const OperatorAbstract* op, bool isLHS, int order);
    static void createIndex(Index* idx, std::vector<int> hierarchy);

    class TensorProduct: public OperatorAbstract{
    private:
        std::vector<UseMatrix*> mats;
        std::vector<int> fillIndicesI;
        std::vector<int> fillIndicesJ;
        std::vector<std::complex<double> >* temporary1;
        std::vector<std::complex<double> >* temporary2;

        void fillIndicesRec(int* fillIndices, int offset, std::vector<int> dims, std::vector<int> filledDims);

    public:
        TensorProduct(std::string Name, const Index* IIndex, const Index* JIndex, std::vector<UseMatrix*> Mats);
        void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

        long applyCount() const;
    };


public:
    const OperatorAbstract* mapToInner;
    const OperatorAbstract* core;
    const OperatorAbstract* mapFromInner;

    static OperatorTucker* truncate(const OperatorAbstract* base);
    void check(const OperatorAbstract* base);
    OperatorTucker(std::string name, std::vector<const OperatorAbstract*> maps);
};


#endif // OPERATOR_TUCKER_H
