// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INITIALSTATE_H
#define INITIALSTATE_H

#include <string>

class OperatorTree;
class Discretization;
class DerivativeFlat;
class Wavefunction;
class Plot;
class Coefficients;
class ReadInput;

class InitialState
{
    static void productFunction(std::string InitialKind, Coefficients & C);
    static void indexOrbital(std::string InitialKind, Coefficients & C);
    static void externalOrbital(std::string InitialKind, Coefficients & C);
    static void subblockEigen(std::string InitialKind, const OperatorTree* Op, Coefficients & C);
public:
    static std::string readKind(ReadInput & Inp);

    /// initial state into wf - switch between various methods
    static Wavefunction get(const Discretization *D, std::string initialKind /** choices: ZERO, Hinitial, atBegin */,
                            int initialN /** N'th excited state */,
                            double tBeg /** time if InitialOper is time-dependent */,
                            const OperatorTree &InitialOper /** initial state can be eigenvector of this operator */,
                            DerivativeFlat * derNew /** time-evolution operator (for initialKind="atBegin") */);
    InitialState();
};

#endif // INITIALSTATE_H
