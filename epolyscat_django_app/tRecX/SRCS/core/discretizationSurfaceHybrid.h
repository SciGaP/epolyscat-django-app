// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONSURFACEHYBRID_H
#define DISCRETIZATIONSURFACEHYBRID_H

#include "operatorAbstract.h"
#include "discretizationSurface.h"

class DiscretizationHybrid;

/// \ingroup Discretizations
/// \brief surface for Neutral (+) Channels discretization
///
/// surface is strictly computed from the Channels part,
/// the neutral does not reach the surface
/// for now only for Neut+Chan type hybrids
class DiscretizationSurfaceHybrid : public DiscretizationSurface
{
public:
    class Map: public OperatorAbstract{
        const OperatorAbstract * surfMap;
    public:
        Map(const OperatorAbstract * SurfMap,const Index* ParentIndex):OperatorAbstract(SurfMap->name,SurfMap->iIndex,ParentIndex),surfMap(SurfMap){}
        void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
        {surfMap->apply(A,*Vec.child(1),B,Y);}
        /// dummy update function for surface
//        void update(double Time, const Coefficients* CurrentVec=0){}
    };

    DiscretizationSurfaceHybrid(const DiscretizationHybrid *Parent, const std::vector<double> &Rad, unsigned int NSurf);
};

#endif // DISCRETIZATIONSURFACEHYBRID_H
