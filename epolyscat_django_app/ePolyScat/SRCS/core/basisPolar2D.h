// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISPOLAR2D_H
#define BASISPOLAR2D_H

#include "basisNdim.h"

class ReadInput;
class BasisPolar2D  : public BasisNdim
{
    std::vector<double> _origin;
    const Index* _idx; // basis in polar coordinates around _origin
    const Index* idxConstruct(double RadMax, int Mmax, int Nrad);
    std::string selectCoor(const std::string Coor) const; ///< select admissible subset
    Index * productIndex(const BasisAbstract *Radial, int Mmax);
protected:
    // all necessary transformations
    double absFactor(const std::vector<double>&CoorNdim) const {return 1.;}
    double nablaFactor(const std::vector<double>&CoorNdim,int I) const {return 0.;}
    Eigen::MatrixXd jacobianToNdim(const std::vector<double>&CoorNdim) const;
public:
    std::vector<double> toCartesian(const std::vector<double>&CoorNdim) const;
    std::vector<double> fromCartesian(const std::vector<double>&CoorNdim) const;
    BasisPolar2D(ReadInput & Inp);
    std::vector<std::complex<double> > operator()(std::vector<double> CoorNdim) const;

};

#endif // BASISPOLAR2D_H
