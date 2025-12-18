// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ODECORSIA_H
#define ODECORSIA_H

#include "odeStep.h"

#include <complex>
#include <vector>
#include <complex>
#include "abort.h"
#include "linSpaceMap.h"
#include "odeSIA.h"


template <class Der,class V>
class OdeCorSIA: public OdeStep<Der,V>{
    OdeSIA<Der,V> *sia;
    V * psi0;
public:
    ~OdeCorSIA(){delete sia,psi0;}
    OdeCorSIA(Der *D,unsigned int MaxKrylov,double EpsSquared)
        :OdeStep<Der,V>("CorSIA",D),
          // NOTE: put somewhat higher accuracy requirements on OdeSIA to ensure correct exponentiations
          sia(new OdeSIA<Der,V>(D,MaxKrylov,std::max(EpsSquared,1.e-12))),
          psi0(new V(D->lhsVector())) {}

    V& step(V &Vec, double Tstart, double Tstep){

        *psi0=Vec;

        OdeStep<Der,V>::derOde->update(Tstart);
        OdeStep<Der,V>::derOde->apply(Tstep/6.,*psi0,0.,Vec);

        sia->step(Vec,Tstart,Tstep);
        OdeStep<Der,V>::nCalls+=sia->currentKrylov();

        sia->step(*psi0,Tstart,Tstep);
        OdeStep<Der,V>::nCalls+=sia->currentKrylov();

        Vec.axpy(1.,*psi0,1.);

        OdeStep<Der,V>::derOde->update(Tstart+Tstep/2.);
        OdeStep<Der,V>::derOde->apply(-Tstep/3.,*psi0,1.,Vec);

        OdeStep<Der,V>::derOde->update(Tstart+Tstep);
        OdeStep<Der,V>::derOde->apply(Tstep/6.,*psi0,1.,Vec);

        OdeStep<Der,V>::nCalls+=3;

        return Vec;
    }

    std::string name() const{return "CorSIA";}
    unsigned int consistencyOrder() const {return 4;}
    double safetyFactor() const {return 0.7;}
    unsigned int currentKrylov(){return sia->currentKrylov();}
    unsigned int maxKrylov(){return sia->maxKrylov();}
};

#endif // ODECORSIA_H
