// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ALGEBRAMULTI_H
#define ALGEBRAMULTI_H

#include <memory>
#include <complex>
#include "algebra.h"

/** \ingroup Functions */

///@brief multi-argument algebra
class AlgebraMulti : public Algebra
{
    mutable std::vector<std::complex<double>> _q;
    size_t _iSingle;
protected:
    std::string _coors; /// set to "*".. for any coordinate, else, e.g., "Eta,Rn", accepts Eta,Rn, but also Eta1,Rn1 etc (numbers must match)
    std::vector<const Algebra*> _arg;
    virtual void checkArgs(const std::vector<std::complex<double> > Q) const;
public:
    static bool isAlgebraMulti(std::string Definition);
    static const AlgebraMulti * factory(const std::string Definition);

    ~AlgebraMulti(){for(auto a: _arg)delete a;}
    AlgebraMulti(std::string Definition /** func(arg1,arg2,...)[:var1,var2..] argN is single-variable Algebra's of varN, varN defaults to QN */
                 ,std::string const Coors);

    /// for use in Algebra tree:
    /// value when all but one coordinate have been set to constants in construction, single variable was renamed to Q
    std::complex<double> val(std::complex<double> Q) const;

    /// value for argument tuple Q
    virtual std::complex<double> valMulti(const std::vector<std::complex<double> > Q) const=0;
};

/// trivial example for a multi-variate function
class AlgebraSum: public AlgebraMulti{
public:
    AlgebraSum(std::string Definition):AlgebraMulti(Definition,"*"){}
    std::complex<double> valMulti(const std::vector<std::complex<double> > Q) const {
        std::complex<double> s(0.);
        for(size_t k=0;k<Q.size();k++)s+=_arg[k]->val(Q[k]);
        return s;
    }
};

/// trivial example for a multi-variate function (asymmetric harmonic oscillator)
class AlgebraExternalHO: public AlgebraMulti{
public:
    AlgebraExternalHO(std::string Definition):AlgebraMulti(Definition,"Eta,Rn"){}
    std::complex<double> valMulti(const std::vector<std::complex<double> > Q) const{
        checkArgs(Q);
        return (1.+0.1*Q[0]*Q[0])*Q[1]*Q[1];
    }
};

/// fill as needed
class AlgebraExternal: public AlgebraMulti{
    std::shared_ptr<const Algebra> trunc;
    double dist,screen;
public:
    AlgebraExternal(std::string Definition);
    std::complex<double> valMulti(const std::vector<std::complex<double> > Q) const;
};

/// fill as needed
class AlgebraCO2Pot: public AlgebraMulti{
    std::shared_ptr<const Algebra> trunc;
    double bCO,screenC,screenO,alfaCharge;
public:
    AlgebraCO2Pot(std::string Definition);
    std::complex<double> valMulti(const std::vector<std::complex<double> > Q) const;
};

/// read rectangular grid from file return linearly interpolated values
class AlgebraNumerical: public AlgebraMulti{
    std::vector<std::vector<double>> _axes;
    std::vector<double> _values;
public:
    AlgebraNumerical(std::string Definition);
    std::complex<double> valMulti(const std::vector<std::complex<double> > Q) const;
};

/// fill as needed
class AlgebraCO2Nuclear: public AlgebraMulti{
    std::shared_ptr<const Algebra> trunc;
    double bCO,screenC,screenO,alfaCharge;
public:
    AlgebraCO2Nuclear(std::string Definition);
    std::complex<double> valMulti(const std::vector<std::complex<double> > Q) const;
};


#endif // ALGEBRAMULTI_H
