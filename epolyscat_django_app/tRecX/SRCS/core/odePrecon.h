// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ODEPRECON_H
#define ODEPRECON_H

#include "odeStep.h"
//#include "derivativeFlat.h"

class OdePrecon : public OdeStep<LinSpaceMap<Coefficients>,Coefficients>
{
    OdeStep<DerivativeFlat,Coefficients>* ode;
public:
    OdePrecon(OdeStep<DerivativeFlat,Coefficients>*Ode):OdeStep(Ode->derivative()),ode(Ode){}
    Coefficients & step(Coefficients &Vec, double Tstart, double Tstep){
         // pre-step action
        ode->derivative()->update(Tstart,std::vector<std::complex<double> >(1,std::complex<double>(Tstart))); // communicate start time
        ABORT("need to implement removeHighEnergies");
        //        ode->derivative()->removeHighEnergies(Vec);
        // make step
        ode->step(Vec,Tstart,Tstep);
        // post-step
//        ode->derivative()->removeHighEnergies(Vec);
        return Vec;
    }
    unsigned int consistencyOrder() const {return ode->consistencyOrder();}
    virtual std::string name() const {return "precon("+ode->name()+")";}
protected:
//    unsigned int safetyFactor(){return ode->safetyFactor();}
};

#endif // ODEPRECON_H
