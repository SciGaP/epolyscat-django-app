// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORRG_H
#define OPERATORRG_H

#include <string>
#include <vector>
#include <complex>

#include "operatorFloor.h"

class UseMatrix;

/** \ingroup OperatorFloors */
/// real general matrix with phase
class OperatorRG: public OperatorFloor
{
    std::vector<double> *dat; ///< pointer to the matrix data
    std::vector<double> *transDat; ///< pointer to transposed storage of the matrix data
    std::vector<double> x,y; ///< temporary reals
    std::shared_ptr<std::vector<std::complex<double>>>phas; ///< pointer to column-wise phases
    std::vector<unsigned int> idx;

protected:
    void axpy(const std::complex<double> &Alfa, const std::complex<double>*X,unsigned int SizX, const std::complex<double>&Beta, std::complex<double> *Y,unsigned int SizY) const;
    void axpyTranspose(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,const std::complex<double> &Beta, std::complex<double>*Y, unsigned int SizY) const;
public:
    OperatorRG(const std::vector<std::complex<double> > &Phas,const std::vector<double> &Dat,unsigned int Rows, unsigned int Cols, std::string Kind);
    void pack(std::vector<int> & Info, std::vector<std::complex<double> >&Buf) const;

    OperatorRG(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf):OperatorFloor("RG"){
        unpackBasic(Info,Buf);
        std::vector<double> rBuf;
        for(unsigned int k=0;k<Info[3];k++){rBuf.push_back(Buf[k].real());rBuf.push_back(Buf[k].imag());}
        dat=addReal(hashString(_rows,_cols),rBuf);
    }
};

#endif // OPERATORRG_H
