// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COORRHOZ_H
#define COORRHOZ_H

#include "coorSystem.h"

/// \ingroup Coordinates
/// \brief Rho.Z, subsystem of cylinder - combines with Eta.R
class CoorRhoZ : public CoorSystem
{
protected:
    std::vector<double> _toRef(const std::vector<double> & RhoZ) const;
    std::vector<double> _fromRef(const std::vector<double> & EtaR) const;
    std::vector<double> _jacRefdCoor(const std::vector<double> & RhoZ) const;
public:
    CoorRhoZ(std::string Name=""):CoorSystem("Rho.Z","Eta.R",Name){}
};

#endif // COORRHOZ_H
