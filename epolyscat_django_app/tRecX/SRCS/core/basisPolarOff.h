// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISPOLAROFF_H
#define BASISPOLAROFF_H

#include "basisNdim.h"
#include "vectorValuedFunction.h"

class ReadInput;

/** \ingroup Basissets */
/// BasisProd placed off-center
class BasisPolarOff : public BasisNdim
{
    std::vector<double> _origin;
    const Index* _idx; // basis in polar coordinates around _origin
    const Index* idxConstruct(double RadMin, double RadMax, int Lmax, int Mmax, int Nrad);
    std::string _strDefinition;
protected:
    // all necessary tranformations (see BasisNdim for explanations)
    std::vector<double> toCartesian(const std::vector<double>&CoorNdim) const;
    double absFactor(const std::vector<double>&CoorNdim) const {return CoorNdim[2];}
    double nablaFactor(const std::vector<double>&CoorNdim,int I) const {if(I==2)return 1.;return 0.;}
    Eigen::MatrixXd jacobianToNdim(const std::vector<double>&CoorNdim) const;

    // for obtaining value
    std::vector<double> fromCartesian(const std::vector<double>&XYZ) const;

public:
    /// construct off-center BasisNdim from input line
    ///
    /// quadGrid will be around off-center origin,
    /// weights, values, and partial derivatives wrt radial functions from 0,0,0, i.e. integration dphi deta dr (w/o r^2 !)
    BasisPolarOff(ReadInput & Inp, std::string Category="PolarOffCenter", int Line=1);
    static Index * productIndex(const BasisAbstract *Radial, int Lmax, int Mmax);
    std::string strDefinition() const {return _strDefinition;}
    std::string str(int Level=0) const {return Level?name():strDefinition();}
    /// true if this overlaps with Other
    bool hasOverlap(const BasisAbstract & Other) const;

    std::vector<std::complex<double> > operator()(std::vector<double> CoorNdim) const;

};

#endif // BASISPOLAROFF_H
