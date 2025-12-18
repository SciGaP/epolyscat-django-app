// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COEFFICIENTSFUNCTION_H
#define COEFFICIENTSFUNCTION_H

#include <string>
#include <complex>
#include <vector>
#include <cfloat>
#include "tree.h"
#include "discretization.h"

class Index;
class Coefficients;


class CoefficientsFunction
{
protected:
    // hold point-wise data
    class Data:public Tree<Data> {
    public:
        Data(){}
        Data(Tree<Data> *Up, Index* I);
        std::vector<std::vector<std::complex<double> > > data;
        Data* at(unsigned int K){return child(K);}
    };
    Data * data;
public:
    CoefficientsFunction():data(0){}
    virtual void multiply(const Coefficients *X,Coefficients *Y,double Time) = 0;
};

class Identity:public CoefficientsFunction{
public:
    void multiply(const Coefficients *X, Coefficients *Y, double Time=DBL_MAX);
};

/// multiply by grid values (if grid, else abort)
class GridValues:public CoefficientsFunction{
public:
    void multiply(const Coefficients *X, Coefficients *Y, double Time=DBL_MAX);
};

// Volkov Phase
class VolkovPhase:public CoefficientsFunction{
public:
    VolkovPhase(const Index* Idx, std::string kLevel);
    VolkovPhase(const Discretization *D, std::string kLevel):VolkovPhase(D->idx(),kLevel){} ///<OBSOLESCENT
    void multiply(const Coefficients *X, Coefficients *Y,double Time);
    static void consistency(std::string Operator); ///< check consistency of operator with Volkov phase
private:
    std::string kLev; // K index on which volkov phase is multiplied
    std::complex<double> Int_iAx,Int_iAy,Int_iAz;
    double t0; // store previous time
    void integrateVectorPotentials(double t);
    void multiplyLocal(const Coefficients *X, Coefficients *Y, Coefficients *phA);
    void setup(Data * Dat, const Index *I, std::vector<std::complex<double> > &Z);

    Coefficients * PhaseAccumulated;
    void UpdateVolkovPhase(Data *Dat, double Time, Coefficients *phA);
};



#endif // COEFFICIENTSFUNCTION_H
