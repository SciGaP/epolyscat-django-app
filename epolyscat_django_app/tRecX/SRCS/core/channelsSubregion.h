// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef UNBOUND_CHANNELS_H
#define UNBOUND_CHANNELS_H

#include <vector>
#include <memory>

#include "coefficientsWriter.h"
#include "str.h"

class DiscretizationSpectralProduct;
class Discretization;
class DiscretizationGrid;
class OperatorTree;
class Coefficients;
class ReadInput;

// #define _CHANNELS_SUBREGION_GRID_

class ChannelsSubregion{
    std::unique_ptr<DiscretizationSpectralProduct> disc;

    // Result before added into average
    std::unique_ptr<Coefficients> coeff;

    // Currently averaged result
    std::unique_ptr<Coefficients> avg;

    // temporary vector
    std::unique_ptr<Coefficients> temp;
    CoefficientsWriter spec;

#ifdef _CHANNELS_SUBREGION_GRID_
    std::unique_ptr<DiscretizationGrid> grid;
    std::unique_ptr<Coefficients> coeffGrid;
    CoefficientsWriter specGrid;
#endif

    static double tStore;

    struct Config{
        std::string hamiltonian;
        std::string kAxis;
        std::vector<std::string> channelAxes;
        std::vector<std::string> partialAxes;

        std::string str(){return Sstr+hamiltonian+kAxis+channelAxes+partialAxes;}
    };
    static Config config;

//    int subregion;
    double tLastAverage;

    // Current average interval. Without meaning if tStore == 0.
    // Initially tBeginAverage == -DBL_MAX
    double tBeginAverage;
    double tEndAverage;

    // Called after vector at Time is placed into coeff
    void write(double Time);

public:
    static void read(ReadInput& Inp);
    ChannelsSubregion(const Discretization* D, const OperatorTree* Ham, std::string Region, ReadInput &Inp);

    void getConfigurations(const Discretization *D, ReadInput &Inp, std::vector<std::string> &configurations);
    static void replaceAll(std::string& s, const std::string &a, const std::string &b);
    bool isCoMCoor(std::string Coordinate);

    void average(const Coefficients* C, double Time);
    void parallelAverage(const Coefficients* C, double Time);
};

#endif // UNBOUND_CHANNELS_H
