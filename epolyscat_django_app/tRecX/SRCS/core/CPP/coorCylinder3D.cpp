// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coorCylinder3D.h"
#include "constants.h"

using namespace std;

std::vector<double> CoorCylinder3D::_fromRef(const std::vector<double> &XYZ) const{

    vector<double> phiRhoZ;
    double rho=sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]);
    if(rho==0.)return {0.,0.,XYZ[2]};

    // note: the two-argument atan takes care of quadrant
    phiRhoZ.push_back(atan2(XYZ[1],XYZ[0]));
    if(phiRhoZ[0]<0)phiRhoZ[0]+=2*math::pi; // map into [0,2pi)
    phiRhoZ.push_back(rho);
    phiRhoZ.push_back(XYZ[2]);
    return phiRhoZ;
}

std::vector<double> CoorCylinder3D::_toRef(const std::vector<double> &PhiRhoZ) const{

    vector<double> xyz;
    xyz.push_back(PhiRhoZ[1]*sin(PhiRhoZ[0]));
    xyz.push_back(PhiRhoZ[1]*cos(PhiRhoZ[0]));
    xyz.push_back(PhiRhoZ[2]);
    return xyz;
}

std::vector<double> CoorCylinder3D::_jacRefdCoor(const std::vector<double> &PhiRhoZ) const{
    double sinPhi=sin(PhiRhoZ[0]);
    double cosPhi=cos(PhiRhoZ[0]);

    vector<double> res;
    res.push_back(-sinPhi*PhiRhoZ[1]);
    res.push_back( cosPhi*PhiRhoZ[1]);
    res.push_back(0.);

    res.push_back(cosPhi);
    res.push_back(sinPhi);
    res.push_back(0.);

    res.push_back(0.);
    res.push_back(0.);
    res.push_back(1.);

    return res;
}
