// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COORETAR_H
#define COORETAR_H

#include "coorSystem.h"

/// \ingroup Coordinates
/// \brief Eta.R  subsystem of polar - combines with Rho.Z
class CoorEtaR: public CoorSystem{
protected:
    std::vector<double> _toRef(const std::vector<double> & Coor) const{ return Coor;}
    std::vector<double> _fromRef(const std::vector<double> & Ref) const{return Ref;}
    std::vector<double> _jacRefdCoor(const std::vector<double> & EtaR) const{return {1.,0.,0.,1.};}
public:
    CoorEtaR(std::string Name=""):CoorSystem("Eta.R","Eta.R",Name){}
};

#endif // COORETAR_H
