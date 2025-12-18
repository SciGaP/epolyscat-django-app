// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INVERSECOULX_H
#define INVERSECOULX_H

#include "inverse.h"
#include "mpiWrapper.h"

class Index;
class OperatorInverseCoulX;

class InverseCoulX:public Inverse
{
    OperatorInverseCoulX *opInvCoulX;
public:
    InverseCoulX(const Index *Idx, unsigned int SubD, unsigned int SuperD, bool BandOvr);
    void apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

    // no correction needed
    void applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{}

    // local versions (not needed)
    void parallelSetup() const{if(MPIwrapper::Size()>1) DEVABORT("cannot be run in parallel");}
    void apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{ABORT("Not implemented!");}
    void applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{ABORT("Not implemented!");}
};

#endif // INVERSECOULX_H
