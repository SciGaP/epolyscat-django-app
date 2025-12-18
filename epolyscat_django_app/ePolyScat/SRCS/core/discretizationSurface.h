// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONSURFACE_H
#define DISCRETIZATIONSURFACE_H

#include <complex>
#include "discretizationDerived.h"
//#include "operator.h"
#include "indexDerived.h"
#include "coefficientsFunction.h"

class OperatorAbstract;

/** \ingroup Discretizations */

/// surface values and derivatives derived from given Discretization
class DiscretizationSurface: public DiscretizationDerived{
    void _construct();
public:
    static std::string prefix; ///< prefix surface output files by this
    static std::string prefixAx; ///< surface axis for coordinate CCC will be prefixAx+CCC
    static std::string readDef(ReadInput& Inp);
    static std::string read(ReadInput& Inp, std::vector<double> &Points, bool &WriteAscii);
    static std::vector<DiscretizationSurface*> all(ReadInput& Inp, std::shared_ptr<Discretization> D,std::vector<double>& surf, std::string Region);

//    static void read(ReadInput & Inp,std::map<std::string,std::vector<double>);
    static DiscretizationSurface * factory(const Discretization *Parent, std::vector<double> &Rad, unsigned int Nsurf, std::string SrcFile="");
    /// construct the full path for the surface file
    static std::string  surfacePath(std::string RunDir, const Index * Idx, std::string Region);


    virtual ~DiscretizationSurface(){}
    DiscretizationSurface(){}
    DiscretizationSurface(const Discretization *Parent, const std::vector<double> & Rad, int NSurf,
                          std::string DerivativeSide="fromInside" /** derivative from side nearer to zero, else boundaryFromBelow, boundaryFromAbove */);
    DiscretizationSurface(const Discretization *Parent, const std::vector<double> &Rad, const std::string SurfaceFile, std::string DerivativeSide="fromInside");

    DiscretizationSurface(std::string SrcFile); // standalon surface - no parent

    std::shared_ptr<const OperatorAbstract> sharedFromParent() const {return _mapFromParent;}

    class IndexS: public IndexDerived {
        // one of the FE axes of the original discretization is converted as follows:
        // - the "continuity level" CCC is replaced by one branch for each surface radius and renamed surfCCC
        //   (it would be good to have the Basis contain the radii, but for now this is just numbered)
        // - the matching floor level is converted to a 2-component Basis ValDer and with value and derivative at the surface
        // - below an exact copy of the remaining index tree is attached
    public:
        IndexS(const Index* I, std::vector<double> Radius, unsigned int NSurf,std::string DerivativeSide="fromInside");
    };
    IndexS * wIdx;

    class MapSurface:public OperatorTree {
    public:
        MapSurface(const Index* SurfI, const Index* FromI);
    };
};

#endif // DISCRETIZATIONSURFACE_H
