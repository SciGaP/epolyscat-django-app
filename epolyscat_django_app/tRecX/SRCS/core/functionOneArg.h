// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef FUNCTIONONEARG_H
#define FUNCTIONONEARG_H
#include <complex>
#include <vector>
#include <map>
#include "tools.h"

/// @brief virtual base class for a range of single-argument functions
///
/// operatorData parses strings and connects to these functions<br>
/// string syntax: functionName[p0,p1,....] where p0,p1,... are parameters
///
/// implements a "Factory" programming pattern
class FunctionOneArg
{
public:
    virtual ~FunctionOneArg(){}
    /// return value for double argument
    virtual std::complex<double> val(double Q) const=0;

    /// check and supplement default parameters
    virtual void defaults()=0;

    /// return value for complex argument
    virtual std::complex<double> val(std::complex<double> Q) const {ABORT(name+" not implemented for complex argument");return 0.;}

    /// return value for "DerivedClassName( double argument)"
    std::complex<double> operator()(double Q) const {return val(Q);}

    /// return value for "DerivedClassName( complex argument)"
    std::complex<double> operator()(std::complex<double> Q) const {return val(Q);}

private:
    /** @name Standard function classes
    */
    /// @cond DEV
    ///@{
    class FunctionPulse;
    class FunctionConst;
    class FunctionSin;
    class FunctionCos;
    class FunctionCos2;
    class FunctionDelta;
    class FunctionGauss;
    class FunctionLegendre;
    class FunctionRGauss;
    ///@}
    /// @endcond
    ///
public:
    static FunctionOneArg* get(std::string Name, bool Abort=true); ///< return pointer to a function defined by name, Abort if not found
protected:
    FunctionOneArg(){}
    FunctionOneArg(std::string Name,std::vector<double>Pars=std::vector<double>(0)):name(Name),pars(Pars){}
    std::vector<double>pars;
    std::string name;
private:
    static std::vector<FunctionOneArg*> avail;
    static std::map<std::string,FunctionOneArg*> list;
    static void add(FunctionOneArg* Func){avail.push_back(Func);} ///< add to list of available functions
    static void standard();                                       ///< set up a list of standard functions
    static void userSupplied();                                   ///< set up a list of user supplied functions
};

#endif // FUNCTIONONEARG_H
