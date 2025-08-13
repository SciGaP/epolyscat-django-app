// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef EIGENTRACE_H
#define EIGENTRACE_H

#include <string>
#include <vector>
#include <complex>
#include <memory>

class Index;
class ReadInput;
class OperatorTree;
class Algebra;
class AsciiFile;

/// trace eigenvalues starting from initial guess for a range of parameters
///
/// (see tutorial/90Floquet)
class EigenTrace
{
    std::vector<std::complex<double>> _eStart;
    std::unique_ptr<Algebra> _parFunc;
    int _nSteps;
    std::string _parUnits;
    double _parMin,_parMax;
    std::unique_ptr<OperatorTree> _op;
    double _ovrErr;
    void trace(std::complex<double> Eguess, std::vector<std::complex<double> > &Par, std::vector<std::complex<double>>& ParEval, AsciiFile &File);
public:
    EigenTrace(ReadInput &Inp);
    bool run(std::string OpDef, const Index *Idx, std::string RunDir);
    void print() const;
};



#endif // EIGENTRACE_H
