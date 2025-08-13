// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OVERLAPDVR_H
#define OVERLAPDVR_H

#include "coefficients.h"
#include "operatorTree.h"

class CoefficientsLocal;
class Index;
/// specialized constructor for DVR overlaps, otherwise strictly OperatorTree
class OverlapDVR : public OperatorTree
{

    void _construct(Coefficients* Diagonal);
    OverlapDVR(Coefficients* Diagonal);
public:
    // note: here we need the Discretization explicitly as no Operator constructor is not optimal
    OverlapDVR(const Index *Idx);
    OverlapDVR(const OperatorTree *OvTree);
};

#endif // OVERLAPDVR_H
