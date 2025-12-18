// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORMAPCHANNELSSURFACE_H
#define OPERATORMAPCHANNELSSURFACE_H

#include "operatorAbstract.h"
#include "odeRK4.h"

class Coefficients;
class Wavefunction;
class Discretization;
class DiscretizationHybrid;
class DiscretizationFactor;
class DiscretizationSpectral;
class DiscretizationSurface;
class DiscretizationGrid;
class TimePropagator;
class RungeKutta4;
class ReadInput;
//class Operator;
class OperatorGradient;
class UseMatrix;
class CoefficientsLocal;

/** \ingroup Structures */
/// map to surface values and derivatives for multiple channels
class OperatorMapChannelsSurface : public OperatorAbstract
{
    const Discretization* _mainDisc;
    const DiscretizationHybrid* hyb;
    DiscretizationSpectral* spec;
    DiscretizationFactor *ion, *ionCh, *freeCh;
    DiscretizationSurface* freeChSurf;
    OperatorGradient * grad;
    Wavefunction *ionWf, *freeWf, *iniState;
    double chanEmax; // maximum ionic channel energy

    TimePropagator* prop;
    OdeStep<LinSpaceMap<CoefficientsLocal>,CoefficientsLocal>*odeLoc;
    OperatorAbstract* chanOp;

    // Initialize to states between Emin and Emax
    void initializeIonWfs();
    double applyThreshold;
    void setProp() const; // set up propagators before first apply

    static void read(ReadInput & In, std::string &ChanOp, std::string &ChanInt,std::string &IonAxes,std::string &Shift, double &Emin, double &Emax);

public:
    static std::string axis(ReadInput & In);
    OperatorMapChannelsSurface(ReadInput& In, const Discretization *Disc);
    void setup(ReadInput & In);
    ~OperatorMapChannelsSurface();

    double energyShift() const;
    void update(double Time, const Coefficients* CurrentVec){_time=Time;if(CurrentVec)return;}
    void apply(std::complex<double> A, const Coefficients &X, std::complex<double> B, Coefficients &Y) const;
    void currentUntwistingMatrix(UseMatrix& res);
    DiscretizationSurface* pointerToChannelSurface();
    void printInfo() const;
};

#endif // OPERATORMAPCHANNELSSURFACE_H
