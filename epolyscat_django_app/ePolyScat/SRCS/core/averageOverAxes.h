// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef AVERAGE_OVER_AXES_H
#define AVERAGE_OVER_AXES_H

#include <vector>
#include <string>

#include "coefficients.h"

class Index;

/**
 * Compute
 * \f[
 * \phi_{km} = \sum_{jl}\psi^*_{jklm}\psi_{jklm}
 * \f]
 * where the indices to be summed over are defined by _axNames.
 *
 * Warning! This class is not extensively tested. Worked for Helium 6D when averaging ober radial indices.
 *
 * This yields the diagonal of the reduced density matrix, so we should have
 *
 * \f[
 * \sum_{km}\phi_{km} = 1
 * \f]
 */
class AverageOverAxes{
    std::vector<bool> averageAxLevel;

    // Suppress warning coefficients construct during timeCritical
    Coefficients temp;

    void setupIndex(const Index* root, Index* parent);
    void averageRec(const Coefficients& src1, const Coefficients& src2, Coefficients& target);
public:
    AverageOverAxes(const Index* idx, std::vector<std::string> _axNames);

    const Index* Idx;
    void average(const Coefficients& src, Coefficients& target);
};


#endif
