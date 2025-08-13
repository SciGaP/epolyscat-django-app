// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ZSLAPACK_H
#define ZSLAPACK_H
#include <complex>
#include "useMatrix.h"

void lapack_zsbgv(UseMatrix &A, UseMatrix &B, UseMatrix & Eval, UseMatrix & Evec,bool Vectors);

class zslapack{
public:
    static void bgv(bool Vectors, int n, int ka, int kb,
                    std::complex<double>* ab, int ldab, std::complex<double>*bb, int ldbb,
                    std::complex<double>* w, std::complex<double>* z, int ldz );
};
#endif // ZSLAPACK_H
