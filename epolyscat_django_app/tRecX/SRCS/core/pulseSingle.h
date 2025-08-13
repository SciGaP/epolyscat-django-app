// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PULSESINGLE_H
#define PULSESINGLE_H

#include <string>
#include "algebra.h"
#include "vectorReal.h"

class PulseSingle{
    friend class Pulse;
public:
    static PulseSingle* factory(std::string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
                                double Ellip=0., double AngleEllip=0., double PolAng=0., double AziAng=0., double DeltaOmega=0.);
    PulseSingle(std::string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
                double Ellip=0., double AngleEllip=0., double PolAng=0., double AziAng=0., double DeltaOmega=0.);
    void resetParameter(std::string Name,double Value);
    std::string str() const;

protected:

    std::string shape;
    double aPeak;     ///< peak envelope value
    double t0;        ///< time of envelope maximum
    double omega;     ///< carrier frequency
    double deltaOmega; ///< linear chirp: om(t) = omega + t * deltaOmega
    double phiCeo;      ///< carrier-envelope offset phase

    double par0,par1; ///< auxiliary parameters (meaning is envelope-dependent)
    double print0,print1; ///< time interval for default printing
    double tBegin,tEnd; ///< beginning and end of pulse (includes infty)

    double Apot(double T);
    double Field(double T);

    double fwhmTransformLimit() const; ///< return chirped FWHM (intensity) of single pulse
    double spectralWidth() const; ///< return chirped FWHM (intensity) of single pulse

    double compMain,compPerp;       ///< relative strength of main and perpendicular components
    VectorReal polarMain,polarPerp; ///< main polarization vector and perpendicular (for elliptic pol)

    double carrier(double T,int Derivative) const;

    virtual double fwhm() const =0; ///< return FWHM (intensity) of single pulse
    virtual double valEnv(double T) const=0;
    virtual double derEnv(double T) const=0;
    virtual double ddEnv(double T) const=0; ///<second derivative of envelope

};

class PulseFlatTop:public PulseSingle {
    Algebra *env,*der;
public:
    ~PulseFlatTop(){delete env; delete der;}
    PulseFlatTop(std::string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
                 double Ellip=0., double AngleEllip=0., double PolAng=0., double AziAng=0., double DeltaOmega=0.);
    double valEnv(double T) const {return env->val(T).real();}
    double derEnv(double T) const {return der->val(T).real();}
    double ddEnv(double T) const {ABORT("not implemented");return 0;}
    double fwhm() const {return par0+par1;}
};

class PulseCosN:public PulseSingle {
    int cosPow;
public:
    PulseCosN(std::string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
              double Ellip=0., double AngleEllip=0., double PolAng=0., double AziAng=0., double DeltaOmega=0.);
    double valEnv(double T) const {
        if(T<-par1 or T>par1)return 0.;
        return std::pow(cos(T*par0),cosPow);
    }
    double derEnv(double T) const {
        if(T<-par1 or T>par1)return 0.;
        return -cosPow*par0*sin(T*par0)*std::pow(cos(T*par0),cosPow-1);
    }
    double ddEnv(double T) const {
        if(T<-par1 or T>par1)return 0.;
        return -cosPow*par0*par0*(std::pow(cos(T*par0),cosPow)
                 -(cosPow-1)*(std::pow(sin(T*par0),2)*std::pow(cos(T*par0),cosPow-2)));
    }
    double fwhm() const;
};

class PulseTrain: public PulseSingle {
    PulseSingle* single;
    double period;
public:
    PulseTrain(std::string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
               double Ellip=0., double AngleEllip=0., double PolAng=0., double AziAng=0., double DeltaOmega=0.);
    double valEnv(double T) const;
    double derEnv(double T) const;
    double ddEnv(double T) const {ABORT("not implemented");return 0;}
    double fwhm() const;
};

class PulseGauss: public PulseSingle {
public:
    PulseGauss(std::string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
               double Ellip=0., double AngleEllip=0., double PolAng=0., double AziAng=0., double DeltaOmega=0.);
    double valEnv(double T) const {return exp(-T*T*par0);}
    double derEnv(double T) const {return -2*T*par0*exp(-T*T*par0);}
    double ddEnv(double T) const {return 0.;}
    double fwhm() const {return sqrt((2*log(2.))/par0);}
};

#endif // PULSESINGLE_H
