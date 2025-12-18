// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisHF.h"
#include "index.h"
#include "coefficients.h"

BasisHF::BasisHF(int NOrbs, const Index* Idx):BasisOrbitalNumerical("*:HF"){
    _orb.resize(NOrbs,Coefficients(Idx));

    for(int i=0;i<NOrbs;i++)_select.push_back(i);
}

void BasisHF::reset(const std::vector<const Coefficients*> & Orbs) const
{
        for(int k=0;k<Orbs.size();k++)_orb[k]=*Orbs[k];
}
