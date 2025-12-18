// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coorRhoZ.h"
#include "constants.h"
#include "str.h"
#include "abort.h"

using namespace std;

std::vector<double> CoorRhoZ::_fromRef(const std::vector<double> &EtaR) const{
    return {EtaR[1]*sqrt(1.-EtaR[0]*EtaR[0]),EtaR[1]*EtaR[0]};
}

std::vector<double> CoorRhoZ::_toRef(const std::vector<double> &RhoZ) const{
    if(RhoZ[1]==0.)return {0.,0.};
    return {RhoZ[2]/RhoZ[1],sqrt(RhoZ[0]*RhoZ[0]+RhoZ[1]*RhoZ[1])};
}

std::vector<double> CoorRhoZ::_jacRefdCoor(const std::vector<double> &RhoZ) const{
    double r=sqrt(RhoZ[0]*RhoZ[0]+RhoZ[1]*RhoZ[1]);
    if(r==0.)DEVABORT(Str("no Jacobian at origin for RhoZ=")+RhoZ);
    double r3=pow(r,3);
    return {-RhoZ[0]*RhoZ[1]/r3,RhoZ[0]/r,RhoZ[0]*RhoZ[0]/r3,RhoZ[1]/r};
}
