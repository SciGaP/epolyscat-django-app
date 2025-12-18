// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INTEGRATETREE_H
#define INTEGRATETREE_H

#include <string>
#include <vector>
#include <map>
#include "coefficients.h"
#include <memory>

class BasisIntegrable;
class Coefficients;

class IntegrationBoundary{
    static std::map<std::string,const IntegrationBoundary*> _list;
public:
    static const IntegrationBoundary* factory(std::string Axis, std::string Definition);
    /// interval(s) to integrate over this coordinate, depending on Q previous coordinates
    virtual const std::vector<std::vector<double> > ranges(const std::vector<double> &Q) const=0;
    /// weight of each integration range (in most cases =1)
    virtual const std::vector<double> rangeWeights(const std::vector<double> &Q) const=0;
};

class IntegrationZoneEta: public IntegrationBoundary{
    std::vector<std::vector<double> > _ranges;
    std::vector<double> _weights;
public:
    IntegrationZoneEta(double LowerTheta,double UpperTheta);
    const std::vector<std::vector<double> > ranges(const std::vector<double> &Q) const {return _ranges;}
    /// where pole is included, take average contributions from either side of the pole
    const std::vector<double> rangeWeights(const std::vector<double> &Q) const
    {return std::vector<double>(_ranges.size(),1./double(_ranges.size()));}
};

/// phi-boundaries of a cone for given phi (see Streak.pdf)
class IntegrationConePhi: public IntegrationBoundary{
    std::vector<std::vector<double> > _ranges;
public:
    IntegrationConePhi(double Phi /** cone axis azimuthal angel */,
                       double Theta /** cone axis polar angel */,
                       double Gamma /** angle from axis to cone boundary */);
    const std::vector<std::vector<double> > ranges(const std::vector<double> &Q) const {return _ranges;}
    const std::vector<double> rangeWeights(const std::vector<double> &Q) const {return {1.,1.};} // all ranges have full weight
};


/// eta-boundaries of a cone for given phi (see Streak.pdf)
class IntegrationConeEta: public IntegrationBoundary{
    const IntegrationBoundary *_phiComp,*_etaComp; // boundaries of the cone's complement (if Gamma>pi/2)
    double theta,gamma,cx,cy,cz,pPart,qPart;
public:
    IntegrationConeEta(double Phi /** cone axis azimuthal angel */,
                       double Theta /** cone axis polar angel */,
                       double Gamma /** angle from axis to cone boundary */);
    const std::vector<std::vector<double> > ranges(const std::vector<double> &Q) const;
    const std::vector<double> rangeWeights(const std::vector<double> &Q) const {return {1.,1.};} // all ranges have full weight
};

/// perform multi-dimensional integral over a Tree-structured function
///
/// (present Coefficients/BasisIntegrable version can be abstracted)
class IntegrateTree
{
    // [constructor] pre-computes coefficients of linear functional on subset of levels
    // integrate(...) forms inner product of weight with integrand where levels match
    // integration ranges are determined from IntegrationBoundary

    std::string _definition;
    std::unique_ptr<Index> _weightIdx,_integralIdx;
    std::unique_ptr<Coefficients>_weights,_integral;
    bool _squareFunctions;

    const BasisIntegrable* lFunc(const Coefficients & C) const; /// return C's function set
    void quadRule(const std::vector<std::vector<double> > &Ranges,const std::vector<double> &RangeWeights, int Pts, std::vector<double> &Quad, std::vector<double> &Weig);
    void quadRule(const std::vector<double> &Range, double RangeWeight, int Pts, std::vector<double> &Quad, std::vector<double> &Weig);
    void setWeights(const Coefficients *Integrand, Coefficients &Weights, std::vector<double> Q);

    // apply linear functional to coefficients
    void integrate(Coefficients &Integral, const Coefficients* Integrand, const Coefficients *Weight);
public:
    IntegrateTree(const std::string Definition, const Coefficients &Func, bool SquareFunctions=false);
    Coefficients integrate(const Coefficients* Integrand); ///< integrate - apply linear functional to Integrand
    const Index* integralIdx() const {return _integralIdx.get();}
};

#endif // INTEGRATETREE_H
