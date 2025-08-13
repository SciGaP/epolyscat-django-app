// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMATDEPENDENT_H
#define BASISMATDEPENDENT_H

#include <string>
#include "qtEigenDense.h"
#include "index.h"
#include "basisMatAbstract.h"


class BasisAbstract;
class UseMatrix;
/// 1-d multiplicative operator that depends on higher axis functions
/// (analoguous to BasisMatMulti, may be abstracted, limited performance and functionality at present)
class BasisMatDependent: public BasisMatAbstract
{
    void construct(std::string Op, const Index * IIndex, const Index* JIndex,std::vector<std::complex<double> > &Diag);
    void valsOnQuadgrid(const std::string TopAx, int Inflate, const Index* Idx,
                        std::vector<double> &Grid, std::vector<double> &Weig, std::vector<std::complex<double> > &Vals);
    std::complex<double> operAtQ(std::complex<double> Q, std::string Op, const Index* IIndex, const Index* JIndex);
    void _construct(std::string Op, const Index * IIndex, const Index* JIndex);
public:
    BasisMatDependent(){}
    BasisMatDependent(std::string Op, const Index * IIndex, const Index* JIndex){_construct(Op,IIndex,JIndex);}
    static std::string modify(const std::string Op, const Index *IIndex, const Index *JIndex); // Operator string name of the depedent axis inserted
};

#endif // BASISMATDEPENDENT_H
