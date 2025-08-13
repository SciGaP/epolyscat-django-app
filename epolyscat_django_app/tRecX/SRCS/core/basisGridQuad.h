// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISGRIDGUAD_H
#define BASISGRIDGUAD_H

#include <vector>

#include "basisAbstract.h"
#include "basisGrid.h"

/// grid with a quadrature rule
class BasisGridQuad : public BasisGrid
{
    static std::map<std::vector<double>,const BasisGridQuad*>_allBasis;
protected:
    BasisGridQuad(const std::vector<double> Mesh,const std::vector<double> Weights):BasisGrid(Mesh),_weights(Weights){}
    BasisGridQuad(const BasisAbstract * Grid);
    std::vector<double> _weights;
public:
    static const BasisGridQuad* factory(const std::string Definition);
    static const BasisGridQuad* factory(const BasisAbstract * Grid);
    static const BasisGridQuad* factory(const std::vector<double> Mesh, const std::vector<double> Weights);
    const BasisGridQuad *append(const BasisAbstract* Grid) const;

    std::string str(int Level=0) const; ///< human-readable description of basis
    std::string strDefinition() const; ///< string fully defines basis
    const std::vector<double> & weights() const {return _weights;}
};

#endif // BASISGRIDGUAD_H
