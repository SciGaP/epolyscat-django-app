#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

#include <memory>
#include "linSpaceVector.h"
#include "timeCritical.h"

#include "str.h"

/// Multivector with Elements of a Hilbert space
template <class Vec>
class MultiVectorHilbert: public LinSpaceHilbert<MultiVectorHilbert<Vec>>{
    std::vector<Vec> c;
    mutable std::shared_ptr<Vec> _intern;
public:
    MultiVectorHilbert(){}

    // for debug only
    MultiVectorHilbert(const MultiVectorHilbert<Vec>& V){
        for(auto cc: V.c){
            c.push_back(cc);
        }
    }

    /// populate with Mult copies of Vec
    MultiVectorHilbert(int Mult, const Vec& V)
    {if(V.idx()==0)DEVABORT("A:0");for(int k=0;k<Mult; k++)c.assign(Mult,V);}
    void assign(int Mult, const Vec& V){if(V.idx()==0)DEVABORT("B:0");for(int k=0;k<Mult; k++)c.assign(Mult,V);}

    // --- LinSpaceHilbert pure virtual functions
    MultiVectorHilbert& axpy(std::complex<double> A, const MultiVectorHilbert & X,std::complex<double> B){
        for(size_t k=0;k<c.size();k++)c[k].axpy(A,X.c[k],B);
        return *this;
    };
    void axpy(std::complex<double> A, const MultiVectorHilbert & X){axpy(A,X,1.);};

    long unsigned int size() const {return (c.size()?c.size()*c[0].size():0);};
    double norm() const {double nrm(0.); for(auto v: c)nrm=std::max(nrm,v.norm());return nrm;};
    std::complex<double> dotProduct(const MultiVectorHilbert & RightHandVector) const{
        std::complex<double> res(0.);
        for(size_t k=0;k<c.size();k++)res+=c[k].dotProduct(RightHandVector.c[k]);
        return res;
    }
    std::complex<double> scalarProduct(const MultiVectorHilbert & RightHandVector) const{
        std::complex<double> res(0.);
        for(size_t k=0;k<c.size();k++)res+=c[k].scalarProduct(RightHandVector.c[k]);
        return res;
    }
    // -------------------------------------

    size_t mult() const {return c.size();};
    Vec& operator()(int I) {return c[I];}
    const Vec& operator()(int I) const {return c[I];}

    /// return C= sum[i] M.c[i]  * A[i]
    Vec & reduceToSingle(std::vector<std::complex<double>> A) const {
        if(A.size()!=c.size())DEVABORT("multivectors do not match size of A");
        if (c.size()==0)_intern.reset(new Vec());
        else if(_intern==0){
            timeCritical::suspend();
            _intern.reset(new Vec(c[0]));
            timeCritical::resume();
        }
        _intern->setToZero();
        for(size_t k=0;k<c.size();k++)_intern->axpy(A[k],c[k]);
        return *_intern;
    }
    /// alternate to above with A vector<double>
    Vec & reduceToSingle(std::vector<double> A) const {
        std::vector<std::complex<double>> cA;
        for(size_t k=0;k<A.size();k++)cA.push_back(A[k]);
        return reduceToSingle(cA);
    }

};

#endif // MULTIVECTOR_H
