// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef COEFFICIENTSMULTI_H
#define COEFFICIENTSMULTI_H

#include <memory>
#include "coefficients.h"
#include "linSpaceVector.h"
#include "timeCritical.h"

// a simple implementation of multiple Coefficients
class OperatorAbstract;
/** \ingroup Coefficients */

/// @brief (incomplete) multiple Coefficients for given Index
class CoefficientsMulti: public LinSpaceHilbert<CoefficientsMulti>
{
    std::vector<Coefficients> c;
    mutable Coefficients _intern;
public:
    CoefficientsMulti(){}

    /// populate with Cols copies of C
    CoefficientsMulti(int Cols, const Coefficients* C);
    CoefficientsMulti(const Index *,unsigned int Cols, std::complex<double> Val=0.);

    CoefficientsMulti& axpy(std::complex<double> A, const CoefficientsMulti & X,std::complex<double> B){
        for(size_t k=0;k<c.size();k++)c[k].axpy(A,X.c[k],B);
        return *this;
    };

    long unsigned int size() const {return (c.size()?c.size()*c[0].size():0);};
    double norm() const {double nrm(0.); for(auto v: c)nrm=std::max(nrm,v.norm());return nrm;};
    CoefficientsMulti& operator*=(std::complex<double> A){for(auto &v: c)v*=A;return *this;};
    std::complex<double> dotProduct(const CoefficientsMulti & RightHandVector) const{
        std::complex<double> res(0.);
        for(size_t k=0;k<c.size();k++)res+=c[k].dotProduct(RightHandVector.c[k]);
        return res;
    }
    std::complex<double> scalarProduct(const CoefficientsMulti & RightHandVector) const{
        std::complex<double> res(0.);
        for(size_t k=0;k<c.size();k++)res+=c[k].scalarProduct(RightHandVector.c[k]);
        return res;
    }

    Coefficients& operator()(int I) {return c[I];}
    const Coefficients& operator()(int I) const {return c[I];}

    /// return C= sum[i] M.c[i]  * A[i]
    Coefficients & reduceToSingle(std::vector<std::complex<double>> A) const {
        if(A.size()!=c.size())DEVABORT("multivectors do not match size of A");
        if (c.size()==0)_intern=Coefficients();
        else if(_intern.size()==0){
            timeCritical::suspend();
            _intern.reset(c[0].idx());
            timeCritical::resume();
        }
        _intern.setToZero();
        for(size_t k=0;k<c.size();k++)_intern.axpy(A[k],c[k]);
        return _intern;
    }

    CoefficientsMulti & apply(std::complex<double> A,const OperatorAbstract & Op, CoefficientsMulti & C);
    UseMatrix innerProduct(const CoefficientsMulti & Rhs) const;
    CoefficientsMulti & operator*=(const UseMatrix & Mat);

    unsigned int vecs() const {return c.size();}
    const Coefficients & vec(unsigned int k) const {return c[k];}
    Coefficients & vec(unsigned int k) {return c[k];}
    void push_back(const Coefficients & C){c.push_back(C);}
    void clear(){c.clear();}
    Coefficients & back(){return c.back();}

    /// orthonormalize, remove singular vectors below threshold EpsRemove
    CoefficientsMulti & orthonormalize(const OperatorAbstract & Overlap,bool Pseudo=false,double EpsRemove=-1.);
};

#endif // COEFFICIENTSMULTI_H
