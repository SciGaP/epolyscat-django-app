// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ODETEST_H
#define ODETEST_H

/// tests for class OdeStep type classes
#include "linSpaceMap.h"
#include "vectorComplex.h"

class VectorComplex;
class TestDer:public VectorComplex,public LinSpaceMap<VectorComplex>{
public:
    TestDer(const VectorComplex & Dat, const VectorComplex &Ht, const VectorComplex & Vec);
    void apply(std::complex<double> A, const VectorComplex &Vec, std::complex<double> B, VectorComplex &Y) const;
    void update(double Time, const std::vector<std::complex<double>> & Parameters=std::vector<std::complex<double>>());
    void update(double Time, const VectorComplex* CurrentVec);
    const VectorComplex & lhsVector() const;
    const VectorComplex & rhsVector() const {return rhsVector();}
private:
    double tCurr;
    VectorComplex dat,vec,ht;
};

class OdeTest{
public:
    static void SIA();
};

#endif // ODETEST_H

