// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __BASICDISC__
#define __BASICDISC__

#include <vector>
#include "discretization.h"
// forward class declarations
class ReadInput;
class UseMatrix;
class Axis;
class IndexConstraint;

/// \ingroup Discretizations
/// \brief standard discretization (product of Axis's)
class BasicDisc: public Discretization {
public:
    BasicDisc(){}
    BasicDisc(std::vector<Axis> Ax, const IndexConstraint* Constraint=nullptr);
    BasicDisc(ReadInput &In, const std::string & Subset="", const IndexConstraint* Constraint=nullptr); ///< take several Axis from input and construct
    BasicDisc(std::string InputFile); ///< create ReadInput from file and construct

    /// create or append to exisint testX.inp and testXY.inp axis inputs for test purposes
    static std::string generateTestInputs(std::string InpFile="", int Dim=1);
    static std::string generateTestInputs(int argc=0, char *argv[]=nullptr, int Dim=1);

    static void Test(); ///< demonstrate general discretization
};

#endif
