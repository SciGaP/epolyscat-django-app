// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORZGSBLOCK_H
#define OPERATORZGSBLOCK_H


#include <string>
#include <vector>
#include <complex>

class UseMatrix;
class Index;
class OperatorZG;

#include "operatorFloor.h"
#include "operatorZG.h"

/** \ingroup OperatorFloors */
/// compress sparse matrix into block (assumes many empty rows or columns)
class OperatorZGSblock: public OperatorFloor
{
    OperatorZG* sub;
    std::vector<unsigned int> colSub,rowSub;
    long applyCount() const { return sub->applyCount(); }
protected:
    void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
              const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const;

    void axpyTranspose(std::complex<double> Alfa, const std::vector<std::complex<double> > & X,
              std::complex<double> Beta, std::vector<std::complex<double> > & Y) const;
public:
    OperatorZGSblock(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf);
    OperatorZGSblock(const UseMatrix *Mat, std::string Kind);
    OperatorZGSblock(const OperatorZGSblock&);

    static unsigned int packCode(){return 2;}
    void pack(std::vector<int>& Info, std::vector<std::complex<double> >&Buf) const;
};

#endif // OPERATORZGSBLOCK_H
