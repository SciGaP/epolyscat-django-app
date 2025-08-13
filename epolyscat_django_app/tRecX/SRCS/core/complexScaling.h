// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COMPLEXSCALING_H
#define COMPLEXSCALING_H
#include <string>
#include <vector>
#include <complex>
#include <climits>
#include <cfloat>
#include <map>
#include <memory>

class ReadInput;
class UseMatrix;

/// \ingroup Coordinates
/// \brief define an impose complex scaling
class ComplexScaling{
    friend class BasisMat;
    friend class BasisSet;
    friend class Axis;
    friend class pmlEta;
    std::string axis;
    std::complex<double> eta;      ///< complex scaling angle
    double _r0upper,_r0lower;        ///< upper/lower scaling radius
    std::string kind;              ///< kinds: ECS, PML (may be extended)
    std::map<std::string,std::shared_ptr<ComplexScaling> >_listScal;
public:
    static std::vector<std::string> names;
    ComplexScaling():axis("ANY"),eta(1.),_r0upper(DBL_MAX),_r0lower(-DBL_MAX),kind("UNSCALED"){} ///< default constructor: no complex scaling
    ComplexScaling(std::string Axis, double theta, double r0upper, double r0lower,std::string Kind="UNSCALED");
    ComplexScaling(ReadInput & in, std::string AxisName); ///< read parameters
    ComplexScaling(std::string Def); ///< from string
    std::string strDefinition() const; ///< complete definition on string
    void coordinates(UseMatrix & z) const;    ///< return scaled coordinates
    std::complex<double> xScaled(double X) const;    ///< return single scaled coordinate point
    std::complex<double> etaX(double X) const; ///< return eta(x)
    std::string str() const; ///< summary string of complex scaling parameters (blank if none)
    double r0up() const {return _r0upper;}
    bool operator==(const ComplexScaling & B) const;
};

#endif // COMPLEXSCALING_H
