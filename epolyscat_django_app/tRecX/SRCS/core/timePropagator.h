// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __TIMEPROPAGATOR__
#define __TIMEPROPAGATOR__
#include <vector>
#include <string>
#include "linSpaceMap.h"
#include "coefficientsLocal.h"
#include "odeRK4.h"
#include "odeSIA.h"

#include "timePropagatorOutput.h"

//forward class declarations
class OperatorAbstract;
class Wavefunction;
class Plot;
class DerivativeFlat;
class ReadInput;

class TimePropagator {

    double h; // actual time step size
    double h_max; // maximal time step size
    double h_fix; // fixed time step (used if !=0)
    double hSum,hSquareSum;
    double desiredPrecision;
    unsigned int successfulTimeSteps;
    unsigned int failedTimeSteps;

    static double _cutE,_applyThreshold; //HACK: this belongs to operators

    TimePropagatorOutput * out;
    OdeStep<LinSpaceMap<Coefficients>,Coefficients> *ode;
    OdeStep<LinSpaceMap<CoefficientsLocal>,CoefficientsLocal> *odeLocal;
    void propagate_intern(Coefficients *C, double &Time, double TEnd); ///< propagate without printing, but possibly writing

public:
    /// constructor for abstract ODE solver and Output
    TimePropagator(OdeStep<LinSpaceMap<Coefficients>,Coefficients> * Ode, TimePropagatorOutput * Out, double Precision, double FixStep);
    TimePropagator(OdeStep<LinSpaceMap<CoefficientsLocal>,CoefficientsLocal> * Ode, TimePropagatorOutput * Out, double Precision, double FixStep);

    /// standard format for reading TimePropagator parameters
    static void read(ReadInput &Inp, double &tBeg, double &tEnd, double &tPrint, double &tStore,
                     double &accuracy, double &cutE, double & FixStep, double &ApplyThreshold,std::string & Method);

    /// propagate and print out
    void propagate(Wavefunction * Wf, double tEnd, const std::string Mode="StartAndStop");

    void fixStep(double StepSize){h_fix=StepSize;}

    std::string str(unsigned int Brief) const;///< string describing the time propagator
    std::string info() const;///< current time-propation info
    void finalize() const; ///< finalize time propagation: flush outputs etc.
    void print(); ///< print definition of time-propagator
};

#endif
