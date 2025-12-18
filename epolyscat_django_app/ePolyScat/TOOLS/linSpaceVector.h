// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef LINSPACEVECTOR_H
#define LINSPACEVECTOR_H

#include "abort.h"
#include <complex>
#include <vector>

/** @defgroup Linalg Linear algebra
 *  @ingroup Tools
 *  \brief implements linear algebra and Hilbert space
 *  @{
 */
/// Element of a linear space
///
/// use by "Curiously Recurring Template Pattern":<br>
/// class MyVector:public LinSpaceVector<MyVector>{...}
template <class V>
class LinSpaceVector
{
public:
    /// this <- A*X + B*this
    virtual V& axpy(std::complex<double> A, const V & X,std::complex<double> B)=0;

    /// vector size (=dimension of linear space)
    virtual long unsigned int size() const=0; // return type long unsingend int conforms with std::vector<...>::size();

    virtual V& operator*=(std::complex<double> A){axpy(0.,*dynamic_cast<V*>(this),A);return *dynamic_cast<V*>(this);};
    virtual void axpy(std::complex<double> A, const V & X){axpy(A,X,1.);};
    virtual V& operator+=(const V & X){return axpy(1.,X, 1.);}
    virtual V& operator-=(const V & X){return axpy(1.,X,-1.);}
};

/// Element from a normed linear space
///
/// use by "Curiously Recurring Template Pattern":<br>
/// class MyVector:public LinSpaceNormed<MyVector>{...}
template <class V>
class LinSpaceNormed:public LinSpaceVector<V>
{
public:
    /// compute and return norm (CAUTION: non necessarily the L2-norm!)
    virtual double norm() const=0;
};

/// Element from a normed linear space
///
/// use by "Curiously Recurring Template Pattern":<br>
/// class MyVector:public LinSpaceHilbert<MyVector>{...}
template <class V>
class LinSpaceHilbert:public LinSpaceNormed<V>
{
public:
    /// dot-product of coefficients (i.e. using complex conjugation for left hand side vectors)
    virtual std::complex<double> dotProduct(const V & RightHandVector) const=0;

    /// scalar product on the Hilbert space (taking into account a metric matrix != 1)
    virtual std::complex<double> scalarProduct(const V & RightHandVector) const =0;

    double normL2sq() const {return std::real(scalarProduct(*dynamic_cast<const V*>(this)));}

};
/** @} */
#endif // LINSPACEVECTOR_H
