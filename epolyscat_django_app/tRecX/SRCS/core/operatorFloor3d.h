// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORFLOOR3D_H
#define OPERATORFLOOR3D_H

#include <string>
#include <vector>
#include <complex>
#include <map>
#include "operatorZD.h"


class ReadInput;
class BasisIntegrable;

typedef std::complex<double> (*potential3d)(const std::vector<std::complex<double> > &);

/** \ingroup OperatorFloors */
/// act on 3-dim basis (OBSOLESCENT)
class OperatorFloor3d : public OperatorZD
{
    // lots of setup, shared by all floor blocks
    //static std::vector<double> OperatorFloor3d::phiGrid;
    //static std::vector<double> OperatorFloor3d::etaGrid;

    static double chargeC;
    static double screenHydrogen;
    static double screenCarbon;
    static double methaneSize;
    static std::string methaneAlign;

    static std::vector<std::vector<std::complex<double> > > potIJ; // potential for IJ-block at quadrature grid
    static std::vector<std::vector<std::complex<double> > > phiBas; // phi-basis at quadrature grid
    static std::vector<std::vector<std::vector<std::complex<double> > > > etaBas; // eta-basis at quadrature grid
    static std::vector<std::vector<std::complex<double> > > invOvr;
    static unsigned int mExpansionSize;
    static unsigned int lExpansionSize;

    static std::complex<double> undefined(const std::vector<std::complex<double> > & Coor){ABORT("define potential by Pot3d: potential");}
    static std::complex<double> identity( const std::vector<std::complex<double> > & Coor){return 1.;}
    static std::complex<double> experimental( const std::vector<std::complex<double> > & Coor){return abs(Coor[2].real())<10.?-1.:0.;}
    static std::complex<double> harmOsc(  const std::vector<std::complex<double> > & Coor){return pow(Coor[2],2)*0.5;}
    static std::complex<double> hydrogen( const std::vector<std::complex<double> > & Coor){return Coor[2]==0.?-DBL_MAX/10.:-1./Coor[2];}
    static std::complex<double> methane(  const std::vector<std::complex<double> > & Coor);
    static std::complex<double> radial( const std::vector<std::complex<double> > & Coor); // general radial potential

    static unsigned int mLev,lLev,femLev; // phi,eta, and FE-DVR and floor levels in hierarchy
    static std::vector<std::vector<std::complex<double> > > basVal(const UseMatrix &Grid, const UseMatrix &Weig, const BasisIntegrable & Bas, unsigned int PotPoints);

    static std::vector<int> quad;
    static std::string potDef;
    static potential3d potOrigin0;
    static std::complex<double> pot3d(const std::vector<std::complex<double> > &);

    static std::vector<double> _potShift;

public:
    OperatorFloor3d(const Index *IIndex, const Index *JIndex);
    static void calculateExpansion(unsigned int lmax);
    static void read(ReadInput & Inp);
    static void print();
    static void setup(const Index * Idx);
};

#endif // OPERATORFLOOR3D_H
