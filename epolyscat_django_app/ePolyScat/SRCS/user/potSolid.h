// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef POTSOLID_H
#define POTSOLID_H

#include <string>
#include "algebra.h"


typedef double(*singleSitePotential)(double Q);

/// collection of functions to compose a solid model
class PotSolid:public Algebra
{
    std::vector<double> lim;
    std::vector<const Algebra*> potSing;
    std::vector<const Algebra*> derSing;
    std::vector<double> lheight,rheight;
    std::vector<double> zero;
    int _der;

    /// return a new PotSolid algebra
    static const Algebra* exportFactory(std::string Term);
    static const Algebra* exportDerivativeFactory(std::string Term);

    /// algebra for single site potential as
    static const Algebra* potSingFactory(std::string Func, int Der=0);
public:

    /// read parameters and register the localFactory at for Algebra::factory
    static void read(ReadInput & Inp){
        PotSolid pot(Inp);
        Algebra::addExternalFactory(exportFactory);
        Algebra::addExternalFactory(exportDerivativeFactory);
    }

    PotSolid(ReadInput &Inp,int Derivative=0);
    std::complex<double> val(const std::complex<double> Q) const;

};

#endif // POTSOLID_H
