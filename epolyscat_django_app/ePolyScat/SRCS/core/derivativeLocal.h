// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DERIVATIVELOCAL_H
#define DERIVATIVELOCAL_H

#include "linSpaceMap.h"
#include "coefficientsLocal.h"
#include "coefficientsGlobal.h"
#include "derivativeFlat.h"
#include "mpiWrapper.h"

class DerivativeLocal : public LinSpaceMap<CoefficientsLocal>
{
    DerivativeFlat * der;
    CoefficientsLocal * model;
    CoefficientsGlobal * globX;
    CoefficientsGlobal * globY;
public:
    DerivativeLocal(){}
    DerivativeLocal(DerivativeFlat* Der):LinSpaceMap<CoefficientsLocal>(Der->applyAlias()),der(Der){
        model=new CoefficientsLocal(Der->idx());
        globX=new CoefficientsGlobal(Der->idx());
        globY=new CoefficientsGlobal(Der->idx());
    }
    void apply(std::complex<double> A, const CoefficientsLocal &X, std::complex<double> B, CoefficientsLocal &Y) const{
           der->apply(A,const_cast<CoefficientsLocal*>(&X),B,Y);
    }
    const CoefficientsLocal & lhsVector() const {return *model;}
    const CoefficientsLocal & rhsVector() const {return *model;}

    void update(double Time, const CoefficientsLocal* CurrentVec=0){ der->update(Time,CurrentVec);}
};

#endif // DERIVATIVELOCAL_H
