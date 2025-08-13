// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef HARMONICS_H
#define HARMONICS_H

#include <string>
#include <vector>

class ReadInput;
class Plot;
class OperatorAbstract;

/// \ingroup Plot
/// \brief compute high harmonic respose from dipole data (uses FFTW)
class Harmonics
{
    int nPoints;
public:
    Harmonics(ReadInput & Inp, std::string Exten, std::string Sep, bool RowWise, double UnitT);
    static std::string dipoleDefinitions(ReadInput & Inp, const std::string &Coor, const std::string & Hamiltonian,
                                         std::vector<std::string> &Names);

    static void compute(int argc, char *argv[]);

    /// find dipole operators in OpList and add to density plots
    static void addDipolesToPlot(std::string DipoleDefinitions, Plot & PlotDef, const std::vector<OperatorAbstract *> &OpList);
};

#endif // HARMONICS_H
