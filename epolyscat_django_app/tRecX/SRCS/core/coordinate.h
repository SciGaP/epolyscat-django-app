// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __COORDINATE__
#define __COORDINATE__
#include <cfloat>
#include <string>
#include <vector>
#include "basisFunction.h"

/// \defgroup Coordinates
/// \ingroup DiscretizationClasses
/// \brief Various single coordinate definitions

/// \ingroup Coordinates

/// \brief list of coordinate systems
class Coordinate{
public:
    static std::map<std::string,Coordinate*> list;
    static void setUp(); ///< create list of available coordinates
    static void cleanUp(); ///< call befor terminating program (to get clean valgrind report)
public:
    enum codes {
        PXi,PEta,
        
        dum,                ///< dummy coordinate
        X,Y,Z,              ///< cartesian \f$(-\infty,\infty)\f$
        Xh,                  ///< cartesian half-axis: \f$[0,\infty)\f$
        Phi,                ///< phi       \f$ [ 0, 2\pi]\f$
        Th,                 ///< theta     \f$ [-\pi/2, \pi/2]\f$
        CTh,                ///< cos(theta)
        R,                  ///< radial direction with square Jacobian
        Rn,                 ///< radial direction without Jacobian
        Rho,                ///< 2d polar radius
        Rhn,                ///< 2d polar radius [srt(rho) absorbed]
        L,                  ///< azimuthal quantum number
        M,                  ///< magnetic quantum number
        Xi,                 ///< prolate spheroidal coordinates: "radial direction" \f$\xi = [1,\infty)\f$
        Eta,                ///< prolate spheroidal coordinates: "azimuthal direction" \f$\eta = [-1,1]\f$
        K,                  ///< prolate spheroidal coordinates: discrete version of eta
        Vec,                ///< "vector", set of indices??
        Idx,                ///< Dummy index coordinate
        Ion,                ///< Ionic index: haCC
        Neut,               ///< Neutral index: haCC
        Orbital             ///< Orbital (expressed in main Discretization)
    };
    /// jacobian for integrals
    enum jacobian {J_one, ///< \f$ \int dq\, f(q) \f$
                   J_val, ///< \f$ \int q\, dq\, f(q) \f$
                   J_squ  ///< \f$ \int q^2\, dq\, f(q) \f$
                  };

    std::string cString; ///< string associated with code
    codes code;           ///<  coordinate code for use in switches etc.
    bool zeroLow,zeroUp; ///< true if lower/upper boundaries =0
    jacobian jaco;       ///< jacobian function code for integrations (turn into function pointer)
    double qmin,qmax;    ///< maximal lower and upper boundary
    std::string function; ///< string name of default function
    bool exactIntegration; ///< default for enforcing exact inegration also in DVR mode
    bool _periodic;       ///< treat as periodic on given (finite) range

public:
    Coordinate():cString("undefined"),code(dum),_periodic(false){}
    /// full definition of coordinate
    Coordinate(std::string String /** name string (same as code) */,
               codes Code         /** code for fast refererence */,
               bool ZeroLow       /** zero boundary at lower end */,
               bool ZeroUp        /** zero boundary at upper end */,
               Coordinate::jacobian Jaco/** Jacobian for integrations */,
               double Qmin        /** lowest possible value */,
               double Qmax        /** highest possible value */,
               std::string Function /** default function kind */,
               bool ExactIntegration /** default einforce exact integration */
               )
        :cString(String),code(Code),zeroLow(ZeroLow),zeroUp(ZeroUp),jaco(Jaco),
          qmin(Qmin),qmax(Qmax),function(Function),exactIntegration(ExactIntegration),_periodic(false){}

    static int automaticOrder(std::string String, int NCoeff); ///< automatically determine order for coordinate axis
    static bool isDiscrete(std::string String); ///< make coordinates contiuous acros axis elements
    static Coordinate fromString(std::string String); ///< return coordinate definition corresponding to string
    static std::string kind(std::string Name);

    Coordinate(std::string String); ///< construct by name
    std::string name() const; ///< return coordinate string
    std::string defaultFunction() const; ///< default choice of fuctions for coordinate
    bool isPeriodic() const {return _periodic;}
    bool isInfinite() const {return qmin<DBL_MAX/2 or qmax > DBL_MAX/2;}
};


#endif // __COORDINATE__
