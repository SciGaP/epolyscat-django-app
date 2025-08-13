// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORBESSELCOULOMB_H
#define OPERATORBESSELCOULOMB_H

#include "operatorFloor.h"

class OperatorAbstract;
class BasisBesselCoulomb;
class OperatorBesselCoulomb:public OperatorFloor
{
    OperatorFloor *floor;
    std::vector<std::complex<double> > bVecI,bVecJ;
    unsigned int mIdxI,mIdxJ;

public:
    OperatorBesselCoulomb(const std::string& TermOper,const BasisBesselCoulomb* IBas,
                          const BasisBesselCoulomb* JBas, std::complex<double> Multiplier);
    void axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
              const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const;
    void pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const {DEVABORT("Not implemented");}
};

#endif // OPERATORBESSELCOULOMB_H
