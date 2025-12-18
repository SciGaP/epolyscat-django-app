// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PULSE_H
#define PULSE_H
#include <complex>
#include <vector>
#include <memory>
#include "vectorReal.h"
#include "algebra.h"
#include "pulseSingle.h"


class ReadInput;
class MultiParam;

/// \ingroup OperatorData
/// \brief laser pulse definitions
class Pulse
        /// a pulse is a vector of 3 components,
        /// each component is the sum of single pulses
{
    friend class FunctionOneArg;
    friend class PulseSingle;
private:
    MultiParam * _parameterRange;

    // store times and matching values to avoid re-computations
    mutable double t0F,t0A,t0D,t0int,intAsqT0;
    mutable VectorReal vecAt0,vecFt0,vecDt0;

    class aSq: public Algebra {
    public:
        aSq(){}
        virtual std::complex<double> val(const std::complex<double> Q) const{return current.Apot2(Q.real());}
    };

    class aField: public Algebra {
        int k;
    public:
        aField(int K):k(K){}
        virtual std::complex<double> val(const std::complex<double> Q) const{return current.Field(Q.real(),k);}
    };


    class aDfield: public Algebra {
        int k;
    public:
        aDfield(int K):k(K){}
        virtual std::complex<double> val(const std::complex<double> Q) const{return current.vecD(Q.real())[k];}
    };

    class aFieldPowN: public Algebra {
        int n;
    public:
        aFieldPowN(int N):n(N){}
        virtual std::complex<double> val(const std::complex<double> Q) const{return std::pow(current.Fabs(Q.real()),n);}
    };

    /// I^N(t) with I(t)...cycle-average of field energy. For multiple frequencies, averager over longest cycle
    class aIntensityPowN: public Algebra {
        const Pulse* _pulse;
        int n;
    public:
        aIntensityPowN(const Pulse* Pulse,int N):_pulse(Pulse),n(N){}
        virtual std::complex<double> val(const std::complex<double> Q) const;
    };

public:
    static Pulse current;
    static void read(ReadInput & Inp,bool Print=false,std::string Kind="Laser"){setCurrent(Inp,Print,Kind);}
    static void setCurrent(ReadInput & Inp,bool Print=false,std::string Kind="Laser");
    static void setCurrent(const Pulse & Pulse,std::string PrintFile="");
    static void printCurrent(std::string File, double T0=0, double T1=0, int Points=200){current.print(File,T0,T1,Points);}

    static double gettBegin(); //{return current.comp[0][0].tBegin;}
    static double gettEnd(); //{return current.comp[0][0].tBegin;}

    static std::complex<double> F0(double Time); ///< electric field
    static std::complex<double> F1(double Time); ///< electric field
    static std::complex<double> F2(double Time); ///< electric field

    static std::complex<double> iAz(double Time); ///< vector potential (time-integral of electric field
    static std::complex<double> iAx(double Time); ///< vector potential (time-integral of electric field
    static std::complex<double> iAy(double Time); ///< vector potential (time-integral of electric field

    static std::complex<double> Axx(double Time); ///< vector potential (time-integral of electric field
    static std::complex<double> Ayy(double Time); ///< vector potential (time-integral of electric field
    static std::complex<double> Azz(double Time); ///< vector potential (time-integral of electric field
    static std::complex<double> Axy(double Time); ///< vector potential (time-integral of electric field
    static std::complex<double> Axz(double Time); ///< vector potential (time-integral of electric field
    static std::complex<double> Ayz(double Time); ///< vector potential (time-integral of electric field



    static std::complex<double> Asqu(double Time); ///< square of vector potential
    static std::complex<double> iAabs(double Time); ///< modulus of vector potential
    static std::complex<double> idThetaA(double Time); ///< vector potential (time-integral of electric field
    static std::complex<double> idThetaF(double Time); ///< vector potential (time-integral of electric field
    static std::complex<double> Fabs(double Time); ///< modulus of field
    static double thetaA(double Time); ///< polar angle of A (for zero field: last at non-zero)
    static double thetaF(double Time); ///< polar angle of A (for zero field: last at non-zero)

    void print(std::string File, double T0=0, double T1=0, int Points=200);
public:
    static bool checkNew; /// debug: globally switch on checking of new pulse
    Pulse():t0F(0),t0A(0),t0D(0),t0int(0),intAsqT0(0){}
    Pulse(ReadInput & Inp, std::string PulseName="Laser");
    void output(std::string Title="", std::string Operator="NONE") const; ///< print pulse parameters
    std::string str(unsigned int Brief=0) const; ///< string representation of pulse parameters
    double beginPrint() const {return print0;} ///< reasonable time for begin print
    double endPrint() const {return print1;}     ///< reasonable time for end print
    double uPonderomotive() const; ///< ponderomotiv potential Amax^2/4
    double intensityMoment(int N) const; ///< pulse energy (N=1) and higher momenta N>1 of intensity: @f$\int dt I^N(t)/2@f$
    double nPhotonCrossSection(double Yield, int N, std::string Units="cm.s") const;
    double apotMax() const; ///< peak vector potential of components
    double omegaMax() const; ///< largest photon energy in pulse
    double omegaMin() const; ///< smallest photon energy in pulse
    double dThetaA (double Time) const; ///< angular velocity for y-axis (d theta / dt)
    double dThetaF (double Time) const; ///< angular velocity for y-axis (d theta / dt)
    double Apot2 (double Time) const; ///< modulus of vector potential

    double Apot (double Time, const unsigned int I) const{return vecA(Time)[I];} ///< i'th component of vector potential
    double Field(double Time,const unsigned int I) const{return vecF(Time)[I];} ///< i'th component of field
    double dField(double Time,const unsigned int I) const{return vecD(Time)[I];} ///< i'th component of field
    void updateParameter(int IComp, std::string ParName, double NewValue);
    void alignFieldZ(double Time); ///< rotate pulse such that field(Time) points in z-direction

    VectorReal vecA(double Time, bool Check=true) const;
    VectorReal vecF(double Time, bool Check=true) const;
    VectorReal vecD(double Time, bool Check=true) const;
    double vecAmax(int Component) const;
    double intAsq(double Time) const;
    double normFsqu(double Time){return vecF(Time,false).normSqu();}
    double normAsqu(double Time){return vecA(Time,false).normSqu();}

    double tBegin() const;
    double tEnd() const;

    // overall information about (multi-)pulse sorted by decreasing wavelength
    std::string polarization(); // z-linear or general (at present)
    double duration();          // approximate overall duration of all pulses
    std::vector<double> wavelengths(); // all wavelengths
    std::vector<double> intensities(); // all inensities

private:
    std::string pulseName;
    std::vector<std::shared_ptr<PulseSingle> >  sing; ///< new single pulses, vectors instead of components (faster)
    double print0,print1; ///< time interval for default printing
    void setFunction(); ///< make three pulse components available as FunctionOneArg
    void parameterRange(ReadInput & Inp); ///< read a parameter range for pulses
    void parameterPts(int ISingle, std::string Name,  std::vector<double*> & Pts); ///< return pointers to the parameter in the three components
};

#endif // PULSE_H
