// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISBESSELCOULOMB_H
#define BASISBESSELCOULOMB_H

#include "basisIntegrable.h"
#include <vector>
#include <complex>

///@brief basis set for the calculation of matrix elements with the Coulomb-Volkov expansion
class BasisBesselCoulomb:public BasisIntegrable
{
    double besselRadius; // use Bessel functions up to this radius
    unsigned int lAngular;
    std::vector<double> kGrid;
    std::vector<double> _invNorms;
    // for the transformation into the standard basis:
    unsigned int _mIdx;
    std::vector<std::complex<double> > bVec;
    BasisIntegrable * _untransformed;
    void valDer(const std::vector<std::complex<double> > &X, const double &K, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der) const;
    void distributeQuadPoints(unsigned int NPoints, int MaxOrder, std::vector<int> &QuadPoints) const;
    void calculateQuadRules(const std::vector<int> &QuadPoints, int &Pos, UseMatrix &X, UseMatrix &W, const double Rmin, const double Rmax) const;
    void initializeTransformation();

    void _construct(double Rc, unsigned int Langular, std::vector<double> KGrid);
public:
    BasisBesselCoulomb(double Rc, double Rx, unsigned int Langular, std::vector<double> KGrid); // standard constructor
    BasisBesselCoulomb(const BasisSetDef &Def); // constructor needed for reading a Bessel-Coulomb basis from an input file

    const BasisIntegrable* pure() const {if(_untransformed==0)return this; return _untransformed;}
    /// values and derivatives at points X, ZeroOutside: set = 0, where X outside range of basis function
    void valDer(const std::vector<std::complex<double> > &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside=false) const;
    bool operator==(const BasisAbstract &other) const;

    std::vector<double> intervals() const {if(_lowBound==besselRadius or _upBound==besselRadius) return {_lowBound,_upBound};return {_lowBound,besselRadius,_upBound};}
    unsigned int mIdx() const {return _mIdx;}
    std::vector<std::complex<double> > bVector() const {return bVec;}
    unsigned int size() const {return kGrid.size();}
    std::string str(int Level=0) const; ///< print basis parameters
    std::string name() const {return "besselCoulomb";}
    bool isGrid() const {return false;}

    void quadRule(int Npoints, std::vector<double> & X, std::vector<double> & W) const;
    unsigned int order()const {return kGrid.size();} ///< order = number of fundamental functions

    void test() const;
};

#endif // BASISBESSELCOULOMB_H
