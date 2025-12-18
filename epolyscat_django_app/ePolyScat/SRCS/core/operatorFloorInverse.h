// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORFLOORINVERSE_H
#define OPERATORFLOORINVERSE_H

#include "operatorFloor.h"
#include "useMatrix.h"

class OperatorFloorInverse:public OperatorFloor
{
    // results from LU factorization; requires "friend class OperatorFloorInverse;" in useMatrix.h!
    UseMatrix *lumat; // L and U matrices
    std::shared_ptr<Eigen::FullPivLU<Eigen::MatrixXcd>> _lu;
    std::vector<int> ipiv; // pivots
    unsigned int subD,superD;
    bool bandOvr; // banded overlap or not

public:
    OperatorFloorInverse(const Index *Idx, unsigned int SubD, unsigned int SuperD, bool BandOvr);
    ~OperatorFloorInverse();
    void axpy(const std::complex<double> & Alfa, const std::complex<double> *X, unsigned int SizX,
                      const std::complex<double> & Beta, std::complex<double> *Y, unsigned int SizY) const;
    void pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const {DEVABORT("Not implemented");}
};

#endif // OPERATORFLOORINVERSE_H
