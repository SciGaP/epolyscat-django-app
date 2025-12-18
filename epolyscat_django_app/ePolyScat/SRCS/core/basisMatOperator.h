// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMATOPERATOR_H
#define BASISMATOPERATOR_H

class BasisAbstract;
class ReadInput;

#include "basisMatMatrix.h"

class BasisMatOperator : public BasisMatMatrix
{
    std::string _operDef;
    std::string _refIndex;
public:
    BasisMatOperator(std::string Kind);
    void setup(const BasisAbstract * IBas, const BasisAbstract * JBas);
};

#endif // BASISMATOPERATOR_H
