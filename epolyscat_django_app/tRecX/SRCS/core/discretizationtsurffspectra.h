// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONTSURFFSPECTRA_H
#define DISCRETIZATIONTSURFFSPECTRA_H

#include "discretizationDerived.h"
#include "index.h"
#include "tree.h"
#include "operatorAbstract.h"
#include "operatorMap.h"

class SurfaceFlux;
class DiscretizationSurface;

/** \ingroup Discretizations */
/// for computing spectra, selected axes converted to momentum grid
class DiscretizationTsurffSpectra: public DiscretizationDerived
{
    void _construct(const Index *IParent, const std::vector<double> Momenta);
public:
    class IndexTsurffSpectra : public Index
    {
    protected:
        std::vector<const Index*> fromIndex; // list of indices from which the present is derived (if emtpy, identical Index tree structures assumed)
    public:
        /// Surface index is converted to momentum index
        IndexTsurffSpectra(const Index * parentI, const std::vector<double> momenta);
    };

    static int _defaultRadialPoints;
    static void computeMomenta(double minEnergy, double maxEnergy, std::vector<double>& momenta, bool KGrid);
    static void readFlags(ReadInput &Inp, int &radialPoints, double &MinEnergy, double &maxEnergy, bool &kGrid, bool AllowFlags);
    static void read(ReadInput &Inp, const Discretization *D, int unboundDOF, std::string &Region, int &NSurf);
    static void checkRegion(int &unbound, std::string &Region, int &NSurf, std::vector<std::string> &infCoorNames, ReadInput &inp);
    static void getNextRegion(int &unbound, std::string &Region, int &NSurf, std::vector<std::string> &infCoorNames, ReadInput &inp);
    static void getRegions(int pos, unsigned int ibeg, int unbound, std::vector<std::string> &infCoorNames, std::vector<std::string> &regions, std::vector<std::string> &names);
    static std::vector<std::string> splitRegion(const std::string s);
    static DiscretizationTsurffSpectra * factoryTsurff(const Index* ParentIdx, ReadInput &Inp, int radialpoints=0);

    /// discretization spectra from surface index
    DiscretizationTsurffSpectra(const Index *FluxIdx, ReadInput &Inp, int RadialPoints=_defaultRadialPoints);
    DiscretizationTsurffSpectra(){}
};

#endif // DISCRETIZATIONTSURFFSPECTRA_H
