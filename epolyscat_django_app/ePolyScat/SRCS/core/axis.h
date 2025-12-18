// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef AXIS_H
#define AXIS_H
#include <climits>
#include <string>
#include <vector>

#include "coordinate.h"
//#include "basisAbstract.h"
#include "complexScaling.h"
#include "basisSetDef.h"

class BasisSetDef;
class BasisAbstract;
class ComplexScaling;

class ReadInput;

/// \ingroup Discretizations
/// @brief A single coordinate axis and the discretization funcitions used on it
class Axis{
    enum AxisKind{FiniteElement};

public:
    static std::string readName(ReadInput & Inp, int Line, std::string Default="BLANK", std::string InputName="Axis");
    static std::string readAlternate(ReadInput & Inp, int Line, std::string InputName,std::string Modify, std::string AlternateModify);

    std::string name;
    Coordinate coor;                 ///< coordinate of the axis
    std::vector<BasisSetDef> basDef; ///< basis function definitions for each element of axis
    std::vector<const BasisAbstract*> bases; ///< alternative to basDef, directly supply bases
    ComplexScaling comsca;           ///< definition of complex scaling (or pml)

    virtual ~Axis(){basDef.clear();bases.clear();}
    Axis(){}
    Axis(const Axis & Ax1, const Axis & Ax2); ///< lossless merge of axes on unscaled, polynomial range
    Axis(std::string Coor, unsigned int Ncoefs, double Qlow, double Qup, std::string Functions="automatic",
         unsigned int Order=1, ComplexScaling Comsca=ComplexScaling()); ///< construct using explicit paramters
    Axis(ReadInput & In, unsigned int Line=1, std::string InputName="Axis", std::string AlternateModify="");  ///< construct single axis from parameters on file
    Axis(const std::string Name,const ComplexScaling Comsca,const BasisSetDef & Def);
    Axis(const std::string Name,const std::vector<const BasisAbstract*> & Bases);

    Axis append(const Axis & App) const; ///< return *this + App
    void appendInPlace(const Axis & App); ///< append App to present axis
    void show(std::string Text="") const; ///< neatly display axis definition
    std::string str(int Brief=0) const; ///< return axis string
    void plot(std::string File, int Points, double QLow=-DBL_MAX, double QUp=DBL_MAX) const; ///< plot basis functions on axis
	
	void remake(bool Deriv);												///< redo axis with specified continuity condition
	void setupXiBasis(std::complex<double> s_k, std::complex<double> q_k);	///< redo axis with specified s_k and q_k for prolate spheroidal coordinates
    void constrain(double Lower, double Upper, unsigned int Order);
    void extendBox(std::vector<double> Box); ///< extend axis to large (unscaled) box, keeping characteristics as closly as possible
    void truncateBox(std::vector<double> Box); ///< remove interval from box (boundaries must coincide with existing element boundaries)
    static void fromFile(ReadInput & In, std::vector<Axis> & ax, std::string Subset="",std::string ReadCategory="Axis"); ///< construct set of axes from input
    static void print(const std::vector<Axis> &ax, std::string File=""); ///< print table with set of axes
    unsigned int maxSize() const; ///< maximun number of coefficients (exact may depend on position in hierarchy)
    double lowerEnd() const; ///< actual lower end of axis (may be \f$ -\infty \f$)
    double upperEnd() const; ///< actual upper end of axis (may be \f$ \infty \f$)
    double lowerRange() const;
    double upperRange() const;
    std::vector<double> elementBoundaries() const;
    double boxsize() const;  ///< Where absorption or laguerre polynomials begin
    unsigned int maxOrder() const; ///< maximal order of all elements on the axis
    unsigned int minOrder() const; ///< minimal order of all elements on the axis
    std::string basString(std::string Kind="brief") const; ///< return human-readable string with basis definitons
    static bool isNameOfAxis(const std::string Name, const std::vector<Axis> Ax, const std::string Mess="");
    bool isElementBoundary(double Val) const;

private:

    void construct(std::string Coor, int Ncoefs,double Qlow, double Qup,std::string Functions,int Order,
                   bool ExactIntegral=false);
    void error(std::string Message);
    void check(); ///< check axis for consistency with definition
};


#endif // AXIS_H
