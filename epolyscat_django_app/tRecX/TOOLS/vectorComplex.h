// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef VECTORCOMPLEX_H
#define VECTORCOMPLEX_H

#include <vector>
#include <complex>
#include "linSpaceVector.h"
#include "tools.h"

/// \ingroup Linalg
/// \brief adds linear space operations to std::vector<complex<double> >
class VectorComplex: public std::vector<std::complex<double> >, LinSpaceHilbert<VectorComplex>
{
public:
    ~VectorComplex(){}
    VectorComplex(){}
    VectorComplex(long unsigned int Dim){resize(Dim);}
    VectorComplex(long unsigned int Dim, std::complex<double> C){assign(Dim,C);}

    VectorComplex& axpy(std::complex<double> A, const VectorComplex & X,std::complex<double> B){
        if(size()!=X.size())ABORT("vector sizes differ");
        const std::complex<double> *x=X.data();
        for(std::complex<double>*y=data();y<data()+size();y++,x++)
            *y=*x*A+*y*B;
        return *this;
    }
    VectorComplex operator*(std::complex<double> A){VectorComplex vA(*this); return vA*=A;}
    VectorComplex & operator*=(std::complex<double> A){for(std::complex<double>*y=data();y<data()+size();y++)*y*=A;return *this;}

    long unsigned int size() const {return std::vector<std::complex<double> >::size();}

    VectorComplex & operator+=(const VectorComplex &X){
        const std::complex<double> *x=X.data();
        for(std::complex<double>*y=data();y<data()+size();y++,x++)*y+=*x;
        return *this;
    }
    VectorComplex & operator-=(const VectorComplex &X){
        const std::complex<double> *x=X.data();
        for(std::complex<double>*y=data();y<data()+size();y++,x++)*y-=*x;
        return *this;
    }

    double norm() const {
        double nrm=0.;
        for(const std::complex<double>*y=data();y<data()+size();y++)nrm=std::max(nrm,std::max(std::abs(y->real()),std::abs(y->imag())));
        return nrm;
    }

    std::complex<double> scalarProduct(const VectorComplex &RightHandVector) const{
        std::complex<double> sum=0.;
        if(size()!=RightHandVector.size())ABORT("vector sizes differ");
        for(const std::complex<double>*y=RightHandVector.data(),*x=data();x<data()+size();y++,x++)sum+=std::conj(*x)*(*y);
        return sum;
    }
    std::complex<double> dotProduct(const VectorComplex &RightHandVector) const {return scalarProduct(RightHandVector);}

    std::string str() const {return tools::str(*this);}
    static void Test(){}
};

#endif // VECTORCOMPLEX_H
