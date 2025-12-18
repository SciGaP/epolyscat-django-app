// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COEFFICENTS_WIRTER_H
#define COEFFICENTS_WIRTER_H

#include <string>

class Coefficients;
class ReadInput;
class AverageOverAxes;
class Wavefunction;

/**
 * Wrapper around Coefficients::write. Given a name and axes to average over(e. g. "Phi.Eta")
 * creates an output folder <name> and <name>/desc.csv listing the values of Index::physical.
 *
 * Can be used from input files:
 *
 * \code
 *      CoefficientsWriter: name,   averageAxes
 *                          coeff,  Rn1.Rn2
 * \endcode
 *
 * Output can be read by SCRIPTS/coefficients_loader.py
 */
class CoefficientsWriter{
private:
    static std::string default_name;
    static std::string default_averageAxes;
    static double default_store;

    /// name == "" <=> disabled
    std::string name;
    std::string averageAxes;
    double store;
    std::string folder;

    bool initial;
    double last_store;
    AverageOverAxes* average;
    Wavefunction* wfOutput;
public:
    static void read(ReadInput& Inp);
    static CoefficientsWriter* instance(); ///< default, i. e. configured instance

    CoefficientsWriter(std::string Name, std::string AverageAxes="", double store=0.);
    void write(double time, const Coefficients& coeffs);

    void disable();
};

#endif
