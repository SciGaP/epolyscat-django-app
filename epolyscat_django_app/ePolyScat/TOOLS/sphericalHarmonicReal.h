// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef SPHERICALHARMONICREAL_H
#define SPHERICALHARMONICREAL_H

#include "abort.h"
#include <complex>
#ifdef _USE_BOOST_
#include "boost/math/special_functions.hpp"
#endif

/// \ingroup Functions
class SphericalHarmonicReal
{
public:
    SphericalHarmonicReal(){}
    double operator()(int L, int M, double Theta, double Phi){
#ifdef _USE_BOOST_
        if(M==0)return std::real(boost::math::spherical_harmonic(L,abs(M), Theta, Phi));
        std::complex<double> plus=boost::math::spherical_harmonic(L, abs(M), Theta, Phi);
//        std::complex<double> minu=boost::math::spherical_harmonic(L,-abs(M), Theta, Phi);
#else
        std::complex<double> plus;
        DEVABORT("no-boost");
#endif
        if(M==0)return real(plus);
        if(M>0)return std::real(plus)*sqrt(2.);
        return std::imag(plus)*sqrt(2.);
    }
};

#endif // SPHERICALHARMONICREAL_H
