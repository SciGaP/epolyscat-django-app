// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
// CODE ADDED BY JOERN STOEHLER

#ifndef BASISDVRSQRT_H
#define BASISDVRSQRT_H

#include <vector>
#include <complex>
#include <memory>

#include "basisIntegrable.h"

#include "basisDvr.h"

class BasisSetDef;
class UseMatrix;
class OrthogonalPolynomial;
/** \ingroup Basissets */

///@brief DVR basis times sqrt(Q)
class BasisDVRSqrt: public BasisDVR
{
public:
    //    static const BasisAbstract* factory(const BasisSetDef & Def);
    BasisDVRSqrt():BasisDVR(){}
    BasisDVRSqrt(const BasisSetDef & Def);
    
    void valDer(const std::vector<std::complex<double> >&X,std::vector<std::complex<double> >&Val,std::vector<std::complex<double> >&Der, bool ZeroOutside) const override;
    std::string str(int Level) const override;
    bool operator==(const BasisAbstract& Other) const override;
};

#endif // BASISDVRSQRT_H
