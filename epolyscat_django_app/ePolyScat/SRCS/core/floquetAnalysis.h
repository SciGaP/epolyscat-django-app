// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef FLOQUETANALYSIS_H
#define FLOQUETANALYSIS_H

#include <string>
#include <vector>

class Index;
class ReadInput;
class OperatorTree;


/// \ingroup Plot
/// \brief read input and set up Floquet analysis
class FloquetAnalysis
{
    std::string select;
    double eMin,eMax,minOvr,maxImag,maxErr;
    std::string dir;
    std::vector<std::string> runs;
    std::string floquetHamiltonian;
public:
    FloquetAnalysis(std::string FloquetHamiltonian, ReadInput & Inp);
    bool run(const Index *Idx);
};

#endif // FLOQUETANALYSIS_H
