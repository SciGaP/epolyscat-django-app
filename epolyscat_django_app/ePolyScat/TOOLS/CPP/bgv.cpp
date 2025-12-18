// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifdef _USE_ZSLAPACK_
#include "zslapack.h"
#include "zslapack_f77.h"
#endif
#include "useMatrix.h"
#include <complex>

using namespace std;

#ifdef _USE_ZSLAPACK_
void zslapack::bgv(bool Vectors, int n, int ka, int kb,
                   std::complex<double> *ab, int ldab, std::complex<double> *bb, int ldbb,
                   std::complex<double> *w, std::complex<double> *z, int ldz)
{
    vector<complex<double > >zwork(3*n);
    int info;
    complex<double>dummy;
    complex<double>*vecs=&dummy;
    char uplow='l';
    char jobz='n';
    if(Vectors){
        jobz='v';
        vecs=z;
    }

    zsbgv_(&jobz,&uplow,&n,&ka,&kb,ab,&ldab,bb,&ldbb,w,vecs,&ldz,zwork.data(),&info);

    if(info!=0){
        cout<<"complex symmetric banded eigensolver (zsbgv) failed, info= "<<info<<endl;
        abort();
    }

}
#endif
