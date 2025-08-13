// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef SCANEIGENVALUE_H
#define SCANEIGENVALUE_H

#include "parameterScan.h"
#include "eigenSubspace.h"

class DerivativeFlat;
//class Operator;
class Coefficients;
class ScanEigenvalue : public ParameterScan
{
    static EigenSubspace * eigSub;
    static std::vector<std::complex<double> > resolv;     // resolvent values
    static std::vector<std::vector<std::complex<double> > > eval;   // eigenvalues for each resolvent value
    static std::vector<std::vector<Coefficients*> > evec;
public:
    ~ScanEigenvalue();
    static void eigenvaluesAtPar(const  std::vector<std::string>& ParName, const std::vector<double> & ParVal, std::vector<double> & Result);
    ScanEigenvalue(ReadInput & Inp);

    void print() const;
    void setEigenSub(const Discretization * D,const OperatorAbstract* H,const OperatorAbstract* H0){eigSub->setup(D,H,H0);}
};

#endif // SCANEIGENVALUE_H
