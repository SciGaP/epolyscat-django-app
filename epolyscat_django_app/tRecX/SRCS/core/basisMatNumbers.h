// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMATNUMBERS_H
#define BASISMATNUMBERS_H

#include "basisMatMatrix.h"

class BasisMatNumbers : public BasisMatMatrix
{
public:
    BasisMatNumbers(const Eigen::MatrixXcd Mat);
    BasisMatNumbers(ReadInput & Inp, int &Line);
};

#endif // BASISMATNUMBERS_H
