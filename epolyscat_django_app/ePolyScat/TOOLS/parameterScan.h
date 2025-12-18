// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARAMETERSCAN_H
#define PARAMETERSCAN_H

#include "readInputRange.h"
#include "asciiFile.h"
/// run through multi-parameter range and store results
class ParameterScan : public MultiParam
{
    typedef void (*ValuesForParameters)(const std::vector<std::string> & Name,const std::vector<double> & Par, std::vector<double > & Result);
    ValuesForParameters valForPar;
    std::string outDir;
public:
    ParameterScan(ReadInput & Inp);
    void setFunction(ValuesForParameters ValForPar){valForPar=ValForPar;}
    void scan();
};

#endif // PARAMETERSCAN_H
