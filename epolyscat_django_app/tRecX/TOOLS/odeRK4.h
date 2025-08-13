// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ODERK4_H
#define ODERK4_H

#include "odeStep.h"

#include <vector>
#include <complex>
#include "abort.h"
#include "unistd.h"
#include "str.h"
#include "parameters.h"

//CAUTION: odeRK4.h is included repeatedly, will cause duplicate timers
//         if multiple time-propagations are run
#include "timer.h"

template<class Der,class V>
class OdeRK4: public OdeStep<Der,V>{
    using OdeStep<Der,V>::derOde; // C++? surprising it does not know where it is derived from
    std::vector<double> a,b,c;
    V *vec0,*kCur,*aux;
public:
    virtual ~OdeRK4(){
//        if(aux!=kCur)delete aux;
//        delete kCur;
//        delete vec0;
    }

    /// classical RK4 scheme
    OdeRK4(Der*D):OdeStep<Der,V>("RK4",D){
        vec0=new V(D->lhsVector());
        kCur=new V(D->lhsVector());

        if(derOde->applyAlias())aux=kCur;
        else aux=new V(D->lhsVector());

        a.assign(4,1./2.);
        b.assign(4,1./3.);
        c.assign(4,1./2.);
        a[0]=0.;a[3]=1.;
        b[0]=1./6.;b[3]=1./6.;
        c[0]=0.;c[3]=1.;
    }

    /// one step of the classical RK4
    V &step(V &Vec, double Tstart, double  Tstep){
        *vec0=Vec;

        for (unsigned int i=0;i<4;i++){
            kCur->axpy(1.,*vec0,Tstep*a[i]);
            derOde->update(Tstart+Tstep*c[i],kCur);      // set time to t0+h*c[i]
            derOde->apply(1.,*kCur,0.,*aux);
            std::swap(kCur,aux);
            Vec.axpy(Tstep*b[i],*kCur,1.);  // y[i+1]=y[i]+h*b[i]*f(t0 + h*c[i+1],y0+h*a[i+1,i]*k[i]), y[4]=final

        }

        OdeStep<Der,V>::nCallsStep=4;
        OdeStep<Der,V>::nCalls+=4;
        return Vec;
    }

    /// return model vector argument of ODEstep
    const V & modelVec() const {return *kCur;}

    unsigned int consistencyOrder() const {return 4;}
    unsigned int nApplyStep() const {return 4;}
    double safetyFactor() const {return 0.1;}
};

#endif // ODERK4_H
