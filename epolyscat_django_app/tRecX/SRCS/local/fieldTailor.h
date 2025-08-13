// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef FIELDTAILOR_H
#define FIELDTAILOR_H
#include <string>
#include <vector>
#include <complex>

#include "linSpaceVector.h"
#include "linSpaceMap.h"
#include "odeStep.h"
#include "vectorReal.h"
#include "readInputRange.h"


#include <dlib/optimization.h>
#include <iostream>

class PhaseSpace: public VectorReal
{
    void setup();
public:
    static std::vector<std::string> compNames;
    PhaseSpace(){setup();resize(6,0.);}
    PhaseSpace(double Qx,double Qy,double Qz,double Px,double Py,double Pz);
    double norm() const{return maxAbsVal();}
    bool inGoing(const std::vector<double> & Poff) const;
    VectorReal q() const {return VectorReal(*this,0,3);}
    VectorReal p() const {return VectorReal(*this,3,6);}
};

typedef dlib::matrix<double,0,1> column_vector;

class AsciiFile;
class ReadInput;
class Pulse;
/// "Tailor" an electric field according to various criteria
class FieldTailor
{    
    std::string kind;// basic selection criterion
    std::string ham; // hamiltonian function for dynamics
    double tMin,tMax; // minimal and maximal times for time-propagation
    double maxStep; // maximal step size in time-propagation
    double tol; // tolerance in propagation
    double rsqApproach; // radious for approach

    bool optimize;
    double rhoBegin;
    dlib::matrix<double>lowLimit;
    dlib::matrix<double>uppLimit;
public:
    ReadInputRange parRange;
private:
    std::string writeComp;
    std::string outDir;
    class Derivative: public LinSpaceMap<PhaseSpace>{
        PhaseSpace lhs,rhs;
    public:
        Derivative(std::string Hamiltonian);
        void apply(std::complex<double> A, const PhaseSpace &X, std::complex<double> B, PhaseSpace &Y) const;
        PhaseSpace & lhsVector() const{return const_cast<Derivative*>(this)->lhs;}
        PhaseSpace & rhsVector() const{return const_cast<Derivative*>(this)->rhs;}
        void update(double Time, const Coefficients* CurrentVec=0);
    };

    class FreeMotionInDipole : public OdeStep<Derivative,PhaseSpace>{
        Pulse* pulse;
        unsigned int _order;
        std::vector<double> quad,weig;
    public:
        FreeMotionInDipole(Pulse* Pul,unsigned int Order);
        PhaseSpace & step(PhaseSpace &Vec, double Tstart, double Tstep);
        unsigned int consistencyOrder() const {return _order;}
    };

    OdeStep<Derivative,PhaseSpace> * dyn; // make a single step by this

    class Trajectory{
    public:
        std::vector<double> time;
        std::vector<PhaseSpace> uVec;
    };

    std::vector<Trajectory> traj;
    Trajectory runSingle(double T0, const  PhaseSpace & UStart) const;
    void zeroIn(Trajectory & Traj) const;
    double startTime() const;
    bool timeExceeded(const Trajectory & Traj) const {return Traj.time.back()>(Traj.time[0]+tMax);}

    VectorReal positionTarget;
    double angleTarget;
    double momentumSquTarget;
    double sigmaMomentumSqu,sigmaAngle;

    void pulseOptimize(std::vector<double> &pars) const;
    double angle(double Time, VectorReal Precol) const;

public:
    ~FieldTailor(){}
    FieldTailor(ReadInput & Inp);
    void run(); // generate a set of trajectories
    Pulse bestMatch() const; // return the pulse with the best match of criteria
    bool stopCriterion(Trajectory &Traj, int & IStop) const;
    bool isRecolliding(const Trajectory & Traj, int BackRange=2) const;

    double operator()(const column_vector& PulseParams) const;

    void print() const;
    void writeData(const AsciiFile &f, std::vector<double> Par, const Trajectory &Traj);
    void writeTrajectories() const;
};

#endif // FIELDTAILOR_H
