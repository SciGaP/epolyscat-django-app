// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INVERSEHYBRID_H
#define INVERSEHYBRID_H
#include<memory>

#include "inverse.h"

class DiscretizationHybrid;
class OperatorVectors;

///@brief Inverse of overlap for hybrid discretization
///
/// H = Ha (+) Hb with non-zero <a|Ovr|b>
///<br> primary use for off-center basis, where overlap remains non-singular
///<br> possible singular vectors are removed from the Ha (smaller) part of the space
///<br> this is mostly not what one wants, if the first part contains a few selected orbitals
///<br> for that case use OperatorSubspace
class InverseHybrid : public Inverse
{
    // initially go through old structure
    std::unique_ptr<OperatorVectors> sInvBA;   // Sb^-1 C
    const OperatorAbstract* sAB;   // C^H
    const OperatorAbstract* sAinv,*sBinv;
    Eigen::MatrixXcd zInv; // Z^-1
    std::unique_ptr<Coefficients> aVec,bVec;

public:
    InverseHybrid(const Index* Idx);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    void apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    void applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

    void apply(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const;
    void apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const;
    void applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const;

    void parallelSetup() const;
};

#endif // INVERSEHYBRID_H
