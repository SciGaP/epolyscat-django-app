// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef BASISMAT_H
#define BASISMAT_H

#include "useMatrix.h"
#include "basisAbstract.h"
#include "abort.h"
#include <map>
#include <cmath>

#include "integrate.h"
#include "algebra.h"

/** @addtogroup OperatorData Operator definitions
 *  @{
*/

/// add a function to the BasisMatFuncs table and check consistency
template<class T> void basisMatFuncSet(T* Func);


/// @brief  abstract base class for user defined functions in matrix calculations
class basisMatFunc {
public:
    virtual ~basisMatFunc(){}

    /// @brief function name as used in operator strings
    ///
    /// example name="One", 1d matrix string "<One[1,2]>" for matrix elements of the characteristic function on interval [1,2]
    std::string name;
    virtual std::complex<double> operator() (std::complex<double>) const
    {ABORT("function undefined for single argument: "+name);} ///< function value
    /// function with vector of arguments
    virtual std::complex<double> operator() (std::vector<std::complex<double> > Q) const{
        if(Q.size()!=1)ABORT("function only defined for single argument: "+name);
        return operator()(Q[0]);
    }

    virtual basisMatFunc & set(const std::string Par){return *this;}         ///< may update parameter
    basisMatFunc(std::string Name):name(Name){basisMatFuncSet<basisMatFunc>(this);}

    virtual bool sameParameters(const basisMatFunc &Other) const {ABORT("multiple construction of userFunction \""+name+"\", need class method sameParameters");}

    /// check whether X has required size (recommended to always use for checking input errors)
    void needSize(std::vector<std::complex<double> > X,unsigned int Size) const {
        if(X.size()!=Size)ABORT(name+" requires argument vector size="+tools::str(Size)+", is="+tools::str(int(X.size())));
    }
};

/// \cond DEV
///< @brief model instantiation of basisMatFunc
class basisMatExp:public basisMatFunc {
    std::complex<double>exponent;
public:
    basisMatExp(std::complex<double> Exponent):basisMatFunc("Exp"),exponent(Exponent){}
    std::complex<double> operator()(std::complex<double> x){return exp(exponent*x);}
    basisMatExp & set(std::string Par=""){if(Par!="")exponent=tools::string_to_complex(Par);return *this;}
};


//! \brief exponential function \f$ \frac{d^n}{(dx)^n}\exp[ikr] \f$ with hindsight to 1d tSURFF
class basisMatExpI: public basisMatFunc {
public:
    basisMatExpI(std::string par = ""): basisMatFunc("ExpI"), der(0), surface(1.) {
        if (par!="") {set(par);}
    }

    virtual std::complex<double> operator() (std::complex<double> k) const {
        return pow(std::complex<double>(0., 1)*k, der)*exp(std::complex<double>(0., surface)*k);
    }

    virtual basisMatFunc& set(const std::string Par) {
        std::vector<std::string> parameters(tools::splitString(Par, ','));
        if (parameters.size()!=2) {
            std::cerr << "Failed to assign new parameters: wrong number of parameters!" << std::endl;
            return *this;
        }
        der=tools::string_to_int(parameters[0]);
        surface=tools::string_to_double(parameters[1]);
        return *this;
    }

    double getSurface() const {return surface;}
    int getDer() const {return der;}
private:
    double surface; // naming with 1d tSURFF in mind
    unsigned int der; ///< order of derivative; 0--value, 1--derivative, -1--integral etc.
};

//! \brief spherical Besselfunction: order of derivative (0--value,1--derivative), order, scaling (abscissa)
class basisMatSpherBessel: public basisMatFunc {
public:
    explicit basisMatSpherBessel(std::string par = ""):
        basisMatFunc("spherBessel"), der(0), order(0), offset(0), surface(1.) {
        if (par!="") { set(par); }
    }
    virtual std::complex<double> operator() (std::complex<double> k) const {
        if (der==0) { // value
            // assume real k for sph_bessel; does not work with complex<double>
            DEVABORT("non-boost");//return sqrt(2./M_PI)*pow(std::complex<double>(0,-1.),order+offset)*boost::math::sph_bessel(order+offset, k.real()*surface);
        }
        else if (der==1) { // derivative; use recurrence relation
            DEVABORT("non-boost");
            return 0.;
            //            return  sqrt(2./M_PI)*pow(std::complex<double>(0,-1.),order+offset)
            //                    *(boost::math::sph_bessel(order+offset, k.real()*surface)/surface*(order+offset)
            //                      -boost::math::sph_bessel(order+offset+1, k.real()*surface)*k.real());
        }
        return 0.;
    }


    virtual basisMatFunc& set(const std::string Par) {
        std::vector<std::string> parameters(tools::splitString(Par, ','));
        if (parameters.size()!=4) {
            std::cerr << "Failed to assign new parameters: wrong number of parameters!" << std::endl;
            return *this;
        }
        unsigned int tempDer =tools::string_to_int(parameters[0]);
        if (tempDer==0 or tempDer==1) {std::swap(der,tempDer);}
        else {std::cerr << "CANNOT HANDLE DERIVATIVE OF ORDER "+parameters[0] << "! DO NOT CHANGE ORDER!\n";}
        order  =tools::string_to_int(parameters[1]);
        offset =tools::string_to_int(parameters[2])/2 +tools::string_to_int(parameters[2])%2;
        surface=tools::string_to_double(parameters[3]);
        return *this;
    }

    unsigned int getDer() const {return der; }
    unsigned int getOrder() const {return order; }
    double getSurface() const {return surface; }
private:
    unsigned int der, order, offset;
    double surface; // with hindsight of using it for tsurff with Volkov states
};
/// \endcond

///// \brief (legacy)
class BasisMat {
public:
    static std::map<std::string,basisMatFunc*> BasisMatFuncs;
};
/// @endcond


template<class T>
void basisMatFuncSet(T* Func){
    std::vector<std::string> akaNames=tools::splitString(Func->name,'|');
    for(unsigned int k=0;k<akaNames.size();k++){
        if(not tools::hasKey(BasisMat::BasisMatFuncs,akaNames[k]))BasisMat::BasisMatFuncs[akaNames[k]]=Func;
        else{
            basisMatFunc* inTable=BasisMat::BasisMatFuncs[akaNames[k]];
            if(inTable->name!=Func->name)return;
            /// note: dynamic_cast may require compiler option on some compilers
            const T * o=dynamic_cast<const T*>(inTable); // need cast, as Other is basisMatFunc
            if(o==0 or not o->sameParameters(*Func))ABORT("inconsistent multiple basisMatFunc="+inTable->name); // cast failed, i.e. classes differ
        }
    }
}
///** @} */ // end group OperatorData


#endif // BASISMAT_H
