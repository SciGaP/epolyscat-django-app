// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COORPOLAR3D_H
#define COORPOLAR3D_H

#include "coorSystem.h"

/// \ingroup Coordinates
/// \brief 3d polar coordinates
class CoorPolar3D: public CoorSystem{
protected:
    std::vector<double> _toRef(const std::vector<double> & PhiEtaR) const;
    std::vector<double> _fromRef(const std::vector<double> & XYZ) const;
    std::vector<double> _jacRefdCoor(const std::vector<double> & PhiEtaR) const;
public:
    CoorPolar3D(std::string Name=""):CoorSystem("Phi.Eta.R","X.Y.Z",Name){}
};

#endif // COORPOLAR3D_H
