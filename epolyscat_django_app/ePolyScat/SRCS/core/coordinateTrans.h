// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COORDINATETRANS_H
#define COORDINATETRANS_H

#include <string>
#include <vector>
#include <map>
#include <complex>


typedef std::vector<double>(*CoordinateMap)(const std::vector<double>& CoorNdim);
typedef std::vector<double>(*JacobianMap)(const std::vector<double>& CoorNdim);
typedef double(*IntegrationFactor)(const std::vector<double>& CoorNdim);
typedef double(*NablaFactor)(const std::vector<double>& CoorNdim,int I);

///\ingroup Coordinates
///@brief collection of coortinate transformation routines (OBSOLESCENT - use CoorSystem and derived classes instead)
class CoordinateTrans
{
public:
    CoordinateTrans();
    unsigned int dim() const; ///< dimension of coordinate system
    std::string & from() const;

    static CoordinateMap toCartesian(std::string To);
    static CoordinateMap fromCartesian(std::string To);
    static JacobianMap jacCartesian(std::string From);
    static IntegrationFactor integrationFactor(std::string Coor);
    static NablaFactor nablaFactor(std::string Coor);

    static double coordinateFrom(std::string New, std::string Old, const std::vector<double>& Vals);

    /// general coordiante transformation (not implemented)

    /// returns vec(r)+vec(shift), where coordinate i/o is in polar form, shift is cartesian
    static std::vector<double> shiftPolar(std::vector<double> PhiEtaR, std::vector<double> CartesianShift);

    static double constOne(const std::vector<double> & Coor){return 1.;}
    static double radius(const std::vector<double> & Coor){return Coor.back();}

    static double nablaOne(const std::vector<double> & Coor,int I){return 0.;}
    static double nablaRadius(const std::vector<double> & Coor,int I){if(I<int(Coor.size())-1)return 0.;return 1.;}

    static std::vector<std::complex<double> > shiftPolar(std::vector<std::complex<double> > PhiEtaR, std::vector<double> CartesianShift);
    static std::vector<double> undefined(const std::vector<double> & In){return In;}
    static std::vector<double> identity(const std::vector<double> & In){return In;}
    static std::vector<double> jacIdentity(const std::vector<double> & In);

    static std::vector<double> fromCylinder3d(const std::vector<double> & In);///< from phi,rho,z to x,y,z
    static std::vector<double> toCylinder3d(const std::vector<double> & Out);///< to phi,rho,z from x,y,z
    static std::vector<double> jacCylinder3d(const std::vector<double> & PhiRhoZ);///< to phi,rho,z from x,y,z

    static std::vector<double> fromPolar3d(const std::vector<double> & In);///< from phi,eta,r to x,y,z
    static std::vector<double> toPolar3d(const std::vector<double> & Out);///< from x,y,z to phi,eta,r
    static std::vector<double> jacPolar3d(const std::vector<double> & PhiEtaR);///< Jacobian d(x,y,z)/d(phi,eta,r) at Phi,Eta,R

    static std::vector<double> fromPolar2d(const std::vector<double> & In);///< from phi,rho to x,y
    static std::vector<double> toPolar2d(const std::vector<double> & Out);///< from x,y to phi,rho
    static std::vector<double> jacPolar2d(const std::vector<double> & PhiRho);///< Jacobian d(x,y)/d(phi,rho) at Phi,Rho
};
#endif // COORDINATETRANS_H
