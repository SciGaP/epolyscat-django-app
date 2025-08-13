// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PLOTSPECTRAL_H
#define PLOTSPECTRAL_H

#include <vector>
#include <string>
#include "plotCoefficients.h"

class DiscretizationSpectral;
class Coefficients;
class OperatorDiagonal;

/// \ingroup Plot
/// \brief plot solution wrt to the spectral representation of a given operator
class PlotSpectral : public PlotCoefficients
{
    const DiscretizationSpectral * spec;
    Coefficients* tmp;
    void intoCols(std::vector<std::complex<double>> Diag, Coefficients * C, std::vector<std::vector<double> > & Cols) const;
public:
    ~PlotSpectral();
    PlotSpectral(const DiscretizationSpectral *Spec);
    /// transform to plot and write to file
    void plot(const Coefficients & C, const std::string & File,
              const std::vector<std::string> & Head=std::vector<std::string>(0),std::string Tag="",bool OverWrite=false) const;
    /// brief name of output file type
    std::string briefName() const {return "sp";}
};

#endif // PLOTSPECTRAL_H
