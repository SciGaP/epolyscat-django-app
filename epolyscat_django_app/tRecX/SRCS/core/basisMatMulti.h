// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMATMULTI_H
#define BASISMATMULTI_H
#include <string>
#include "qtEigenDense.h"
#include "index.h"
#include "basisMatAbstract.h"

class BasisAbstract;
class UseMatrix;
/// multiplicative operator (only in DVR, diagonal)
class BasisMatMulti: public BasisMatAbstract
{
    void construct(std::string Op, const Index * IIndex, const Index* JIndex,std::vector<std::complex<double> > &Diag);
protected:
    void _construct(std::string Op, const Index *IIndex, const Index *JIndex);
public:
    BasisMatMulti(){}
    BasisMatMulti(std::string Op, const Index * IIndex, const Index* JIndex){_construct(Op,IIndex,JIndex);}
};

#endif // BASISMATMULTI_H
