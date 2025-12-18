// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ODEeuler_H
#define ODEeuler_H

#include "odeStep.h"

#include <vector>
#include <complex>
#include "abort.h"
#include "unistd.h"
#include "str.h"
#include <memory>

/// explicit Euler
template<class Der,class V>
class OdeEuler: public OdeStep<Der,V>{
    using OdeStep<Der,V>::derOde; // C++? surprising it does not know where it is derived from
    std::unique_ptr<V> lhs;
public:
    virtual ~OdeEuler(){}

    /// classical euler scheme
    OdeEuler(Der*D):OdeStep<Der,V>("euler",D){lhs.reset(new V(derOde->lhsVector()));}

    /// one step explicit euler
    V &step(V &Vec, double Tstart, double  Tstep){
        derOde->update(Tstart);
        *lhs=Vec;
        derOde->apply(Tstep,*lhs,1.,Vec);
        return Vec;
    }
    /// return model vector argument of ODEstep
    const V & modelVec() const {return*derOde->rhsVector();}

    unsigned int consistencyOrder() const {return 1;}
    unsigned int nApplyStep() const {return 1;}
    double safetyFactor() const {return 0.5;}
};

#endif // ODEeuler_H
