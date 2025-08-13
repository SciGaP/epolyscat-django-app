// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef VISUALIZE_ANGULAR_DISTRIBUTION
#define VISUALIZE_ANGULAR_DISTRIBUTION

#include <map>
#include <string>

class Index;
class Coefficients;

class VisualizeAngularDistribution{
    std::map<const Index*, double> data;

    std::string x_axis;
    std::string y_axis;

public:
    VisualizeAngularDistribution(const Index* Root);
    VisualizeAngularDistribution(const Coefficients* Coeff);

    VisualizeAngularDistribution& withAxes(std::string XAxis, std::string YAxis){
        x_axis = XAxis;
        y_axis = YAxis;
        return *this;
    }

    void print();

    static void printAllPossible(const Index* Idx);
};


#endif
