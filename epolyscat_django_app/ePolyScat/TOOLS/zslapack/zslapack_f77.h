// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ZSLAPACK_F77_H
#define ZSLAPACK_F77_H
#include "stdlib.h"
#include <complex>

/// direct interfaces to the Fortran77 zslapack routines
/// do not use unless you know what you are doing
/// standard interfaces are in "zslapack.h"
extern "C" {

void zsbgv_(const char *vect,const char *uplow,
            const int *n,const int *ka,const int *kb,
            std::complex<double> *A, const int *lda,
            std::complex<double> *B, const int *ldb,
            std::complex<double> *Eval,
            std::complex<double> *Evec, const int *ldv,
            std::complex<double> *zwork,int *info);
}

#endif // ZSLAPACK_F77_H
