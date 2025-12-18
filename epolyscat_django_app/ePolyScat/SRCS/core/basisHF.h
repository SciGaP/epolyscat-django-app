// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISHF_H
#define BASISHF_H

#include "basisOrbitalNumerical.h"

class BasisHF : public BasisOrbitalNumerical
{
public:
    BasisHF(int NOrbs, const Index* Idx);//:BasisOrbital("HF"){_orb.resize(NOrbs);}
    void generateOrbitals(const Index* Idx=0){} ///< here dummy, construction directly
    void reset(const std::vector<const Coefficients*> &Orbs) const;
};

#endif // BASISHF_H
