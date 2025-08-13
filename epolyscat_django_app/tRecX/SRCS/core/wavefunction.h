// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __WAVEFUNCTION__
#define __WAVEFUNCTION__

#include "tools.h"

//forward class declarations
class Coefficients;
class Discretization;
class Index;


/** \ingroup Coefficients */

/// \brief Coefficients plus a time-parameter (obsolescent)
class Wavefunction{
    bool _map; // admittedly not best practice, but historical...
public:

    Wavefunction():time(0.),coefs(0),_map(false){}
    Wavefunction(const Wavefunction& other);
    Wavefunction(const Discretization* D, double Time=0.);
    Wavefunction(const Index* Idx, double Time=0.);
    Wavefunction(double Time, Coefficients *C);

    static const Wavefunction map(double Time, const Coefficients* C);

    virtual ~Wavefunction();

    //natural operations
    Wavefunction & operator=(const Wavefunction &rhs);  //NOTE: DOES set time equal !
    Wavefunction & operator*=(std::complex<double> x);  //NOTE: does not set time equal !
    Wavefunction & operator+=(const Wavefunction &rhs); //NOTE: does not set time equal !
    Wavefunction & operator-=(const Wavefunction &rhs); //NOTE: does not set time equal !
    void setToZero();
    void setToFunction(std::string Function); ///< initialize with function
    void setToConstant(std::complex<double> c);
    void setToRandom(int seed = 1337); ///sets coefficients to random complex numbers with a given seed
    void makeContinuous(); /// imposed continuity conditions
    double maxCoeff();
    void show(int Digits=0);

    //the central obejcts
    double time;
    Coefficients* coefs;

    static double readLastTime(std::ifstream &Inp);
    void  read(std::ifstream & stream,bool header); ///< binary read
    void write(std::ofstream & stream,bool header) const; ///< binary write
    void print(std::ofstream & stream,std::string Header) const; ///< formatted write

private:

};

#endif
