// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ODEARNRK_H
#define ODEARNRK_H

#include "odeStep.h"

#include <complex>
#include <vector>
#include <complex>
#include "abort.h"
#include "linSpaceMap.h"
#include "arnoldi.h"
#include "odeSIA.h"
#include "odeRK4.h"

template <class Der, class V>
class DerPrecon: public LinSpaceMap<V>
{
    OdeSIA<Der,V> *sia;
    Der * der;
    V* phi;
    double tCurr,tCenter;
public:
    ~DerPrecon(){delete sia,phi;}
    DerPrecon(OdeSIA <Der,V> * Sia,Der *D)
        :sia(Sia),der(D),phi(new V(D->lhsVector())){}

    void apply(std::complex<double> A, const V &Vec, std::complex<double> B, V &Y) const
    {
        // U0-propagate to current time
        *phi=Vec;
        der->update(tCenter);
        sia->step(*phi,tCenter,tCurr-tCenter);

        // apply H(t)-H(t0)
        der->update(tCurr);
        der->apply(1.,*phi,0.,Y);
        der->update(tCenter);
        der->apply(-1.,*phi,1.,Y);

        // back-propagate from current time
        sia->step(Y,tCenter,tCenter-tCurr);
    }

    void update(double Time,const std::vector<std::complex<double> > & Par=std::vector<std::complex<double> >()) {
        tCurr=Time;
        if(Par.size()==1)return;
        tCenter=Par[0].real();
    }

    const V & lhsVector() const {return *phi;}
    const V & rhsVector() const {return *phi;}

};

template <class Der,class V>
class OdeArnRK: public OdeStep<Der,V>{
    DerPrecon<Der,V> * derPrecon;
    OdeRK4<DerPrecon<Der,V>,V> *rk4;
    OdeSIA<Der,V> *sia;
public:
    ~OdeArnRK(){delete rk4,sia,derPrecon;}
    OdeArnRK(Der *D,unsigned int MaxKrylov,double EpsSquared):OdeStep<Der,V>("ArnRK",D)
    {
        sia=new OdeSIA<Der,V>(D,MaxKrylov,EpsSquared);
        derPrecon=new DerPrecon<Der,V>(sia,D);
        rk4=new OdeRK4<DerPrecon<Der,V>,V>(derPrecon);
    }

    V& step(V &Vec, double Tstart, double Tstep){
        std::vector<std::complex<double> > par;
        par.push_back(Tstart+Tstep/2.);
        derPrecon->update(Tstart,par);
        rk4->step(Vec,Tstart,Tstep);
        sia->step(Vec,Tstart,Tstep);
       return Vec;
    }

    std::string name() const{return "ArnRK";}
    unsigned int consistencyOrder() const {return 4;}
    double safetyFactor() const {return 0.5;}
};

#endif // ODERK4_H
