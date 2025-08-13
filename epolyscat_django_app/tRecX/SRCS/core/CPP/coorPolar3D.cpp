// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coorPolar3D.h"
#include "constants.h"
#include "str.h"

using namespace std;

std::vector<double> CoorPolar3D::_fromRef(const std::vector<double> &XYZ) const{

    vector<double> phiEtaR;
    double r=sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2]);
    if(r==0.)return vector<double>(3,0.);

    // note: the two-argument atan takes care of quadrant
    phiEtaR.push_back(atan2(XYZ[1],XYZ[0]));
    if(phiEtaR[0]<0)phiEtaR[0]+=2*math::pi;
    phiEtaR.push_back(XYZ[2]/r);
    phiEtaR.push_back(r);
    return phiEtaR;
}

std::vector<double> CoorPolar3D::_toRef(const std::vector<double> &PhiEtaR) const{

    vector<double> xyz;
    double sinThe=sqrt(1.-PhiEtaR[1]*PhiEtaR[1]);
    xyz.push_back(PhiEtaR[2]*sinThe*cos(PhiEtaR[0]));
    xyz.push_back(PhiEtaR[2]*sinThe*sin(PhiEtaR[0]));
    xyz.push_back(PhiEtaR[2]*PhiEtaR[1]);
    return xyz;
}

std::vector<double>  CoorPolar3D::_jacRefdCoor(const std::vector<double> &PhiEtaR) const{

    /// ( dx/dphi  dx/deta  dx/dr )<br>
    /// ( dy/dphi  dy/deta  dy/dr )<br>
    /// ( dz/dphi  dz/deta  dz/dr )

    double sinPhi=sin(PhiEtaR[0]);
    double cosPhi=cos(PhiEtaR[0]);
    double cosThe=PhiEtaR[1];
    double sinThe=sqrt(1.-PhiEtaR[1]*PhiEtaR[1]);

    vector<double> res;
    res.push_back(-sinPhi*sinThe*PhiEtaR[2]);
    res.push_back( cosPhi*sinThe*PhiEtaR[2]);
    res.push_back(0.);

    res.push_back(-cosPhi*cosThe/sinThe*PhiEtaR[2]);
    res.push_back(-sinPhi*cosThe/sinThe*PhiEtaR[2]);
    res.push_back(PhiEtaR[2]);

    res.push_back(cosPhi*sinThe);
    res.push_back(sinPhi*sinThe);
    res.push_back(cosThe);

    return res;
}
