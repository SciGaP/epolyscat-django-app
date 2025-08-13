// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORFPRODUCT_H
#define OPERATORFPRODUCT_H

#include <vector>
#include <complex>
#include <operatorFloor.h>


/// product floor: Y  = prod[i] Floors[i] X
class OperatorFProduct : public OperatorFloor
{
    std::complex<double> _factor;
    std::vector<const OperatorFloor*> _floors;
    bool _view;
    mutable std::vector<std::vector<std::complex<double>>> _tmps;

public:
    ~OperatorFProduct();
    OperatorFProduct(std::complex<double> Factor,std::vector<const OperatorFloor*> Floors, bool View=false);
    void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
              const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY)const;
    void pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const ;
};

#endif // OPERATORFPRODUCT_H
