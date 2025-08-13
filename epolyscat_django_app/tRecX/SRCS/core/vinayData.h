// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef VINAYDATA_H
#define VINAYDATA_H

#include "tools.h"
#include <string>
#include <vector>

class Coefficients;
class Index;
class BasisAbstract;

///@brief reads spectral amplitudes from a ascii file
///
class VinayData
{
    /// structure of the Ascii ampl_x files
    ///
    /// - expected in subdirectory "polarAngle_phiCEO" such as 45_pi2
    ///   where first part is degrees, second part is for pi/2 etc.
    ///
    /// - kGrid does not need to be exact final plot grid, there is an option
    ///   to interpolate and/or sample from a larger grid (introduced as initially
    ///   we had k-Grids that were way to fine over a very large k-range
    ///
    /// - doubles should be given as sufficiently many (e.g. 14) ascii digits
    ///
    /// line format   name         meaning
    ///<br> 1:   1*int     code          signals further format, at present code=1 expected
    ///<br> 2:   3*int     nPhi nEta nK  number of grid points in respective direction
    ///<br> 3: nPhi*double Phi-grid, should be equidistant 2*pi*n/nPhi, n=0,...,nPhi-1
    ///<br> 4: nEta*double Eta-grid,
    ///<br> 5: nEta*double Eta-weights
    ///<br> 6: nK*double   k-Grid
    ///<br> 7: nPhi*nEta*nK complex ... blank-separated amplitude values,
    ///          (2.09655384100621e-06,1.21639224558644e-06) (...,...) (...) ...
    static std::vector<const Index*> _allIdx;
    Coefficients * _c;
    Index* _idx;

    /// convert a Phi.Eta.kRn quadrature grid into an Index
    Index* getIndex(std::vector<std::vector<std::vector<double> > > Grids, std::vector<const BasisAbstract*> Bas, int Level);
    /// read single Ascii line with blank-separated floats into vector<double>
    void readVec(std::ifstream & Inp, int Size, std::vector<double> & Vec);
    double _phiCEO;
    double _polarAngle;
    std::string _channel;
    std::string _fullGrid;
    std::string _sampleGrid;
    std::vector<double> _sideBand;
public:
    ~VinayData();

    ///@brief FileName encodes alignment and phiCEO
    ///
    /// expect format VinayData/45_pi4 for a calculation at 45 degrees alignement and carrier-envelope offset of pi/4
    /// 45_pi4/ampl_x, x=0,...,5 has amplitudes for the individual channels
    VinayData(std::string &FileName);
    const Coefficients* ampl() const {return _c;}
    const Index* idx() const {return _idx;}
    double phiCEO() const {return _phiCEO;}
    double polarAngle() const {return _polarAngle;}
    std::string channel(){return _channel;}
    std::string fullGrid(){return _fullGrid;}
    std::string sampleGrid(){return _sampleGrid;}
    std::vector<double> kSideBand() const {return _sideBand;}
};

#endif // VINAYDATA_H
