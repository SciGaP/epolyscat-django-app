// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COORCYLINDER3D_H
#define COORCYLINDER3D_H

#include "coorSystem.h"

/// \ingroup Coordinates
/// \brief cylinder coordinates Phi.Rho.Z
class CoorCylinder3D : public CoorSystem
{
protected:
    std::vector<double> _toRef(const std::vector<double> & PhiRhoZ) const;
    std::vector<double> _fromRef(const std::vector<double> & XYZ) const;
    std::vector<double> _jacRefdCoor(const std::vector<double> & PhiRhoZ) const;
public:
    CoorCylinder3D(std::string Name=""):CoorSystem("Phi.Rho.Z","X.Y.Z",Name){}
};

#endif // COORCYLINDER3D_H
