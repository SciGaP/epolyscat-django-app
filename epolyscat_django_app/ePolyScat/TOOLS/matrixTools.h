// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef MATRIXTOOLS_H
#define MATRIXTOOLS_H

#include "toolsHeader.h"

namespace tools{

/// UseMatrix style string
std::string strMatrix(std::complex<double>* Mat,int Rows,int Cols,int Digits=3,std::string Text="");

/// block-wise string: largest in magnitude element
std::string strMatrixBlock(std::complex<double>* Mat,int Rows,int Cols,
                           int RowBlock=1 /** fixed block size, (except last) */,
                           int ColBlock=0,/** defaults to ColBlock=RowBlock */
                           int Digits=3,std::string Text="");

/// block-wise string: show largest in magnitude element
std::string strMatrixBlock(std::complex<double>* Mat,
                           const std::vector<int>IBlock /** row block sizes */,
                           const std::vector<int>JBlock /** column block sizes */,
                           int Digits=3,
                           std::string Text="");

/// set small real or imaginary part = 0 (smallness relative to column and row magnitude)
void purgeMatrix(std::complex<double>* Mat, int Rows, int Cols, double Eps=1.e-12, double EpsAbs=1.e-14);

}
#endif // MATRIXTOOLS_H
