// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ODESTEP_H
#define ODESTEP_H

#include <vector>
#include <complex>
#include "abort.h"

#include "timer.h"
TIMER(step0,)
TIMER(step1,)
TIMER(step2,)

/// ODE solver template class
///
/// Der ...derivative operator - requirements as LinSpaceMap<br>
/// V .....vector - requirements as LinSpaceNormed
template<class Der ,class V >
class OdeStep
{
protected:
    /// a safety factor for computing new step size (default provided)
    virtual double safetyFactor() const {return 0.6;}

    /// derivative operator - must comply with abstract base class LinSpaceMap
    Der * derOde;
    unsigned int nCalls;
    unsigned int nCallsStep;
    std::string _name;
    unsigned int _consistency;
public:

    virtual ~OdeStep(){}
    OdeStep(std::string Name, Der* D):_name(Name),derOde(D),nCalls(0),nCallsStep(0){}
    /// single step of an ODE solver: overwrite Vec with its value at time Tstart+Tstep, return reference to Vec
    virtual V & step(V & Vec, double Tstart, double Tstep)=0;

    /// return consistency order of method, abort if undefined
    virtual unsigned int consistencyOrder() const {ABORT("no consistency order defined (needed for step-size prediction)");return 0;}

    /// access to derivative that is being used
    virtual Der* derivative(){return derOde;}

    /// return error vector obtained from comparing Tstep with 2 x Tstep/2
    virtual V & stepError(V & Vec, double Tstart, double Tstep, V & Err){


        Err=Vec;
        step(Err,Tstart,Tstep);
        step(Vec,Tstart,          Tstep*0.5);
        step(Vec,Tstart+Tstep*0.5,Tstep*0.5);

        Err-=Vec;
        return Err;
    }

    /// advance by two half-steps and return error vector obtained from comparing Tstep with 2 x Tstep/2
    virtual V & stepError(std::vector<V> & Vec, double Tstart, double Tstep, V & Err){
        Err=Vec.back();
        step(Err,Tstart,Tstep);

        Vec.push_back(Vec.back());
        step(Vec.back(),Tstart,          Tstep*0.5);
        Vec.push_back(Vec.back());
        step(Vec.back(),Tstart+Tstep*0.5,Tstep*0.5);
        Err-=Vec.back();
        return Err;
    }

    /// get new step size and decide whether to accept Error.norm()
    bool acceptStep(double Epsilon, double StepCurrent, double & StepNext, double Error) const {
        if(Error>0.)
            StepNext=StepCurrent*std::pow(this->safetyFactor()*Epsilon/Error,1./double(consistencyOrder()+1));
        else
            StepNext=StepCurrent*2;
        return Error<Epsilon;
    }

    std::string name() const {return _name;}

    unsigned int nApplyFunction() const {return nCalls;}
    virtual unsigned int nApplyStep() const {return nCallsStep;}
    virtual std::string info() const {return name()+": calls="+std::to_string(nCalls);}

};

#endif // ODESTEP_H
