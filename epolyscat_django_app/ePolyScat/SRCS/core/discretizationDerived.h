// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __DISCRETIZATION_DERIVED__
#define __DISCRETIZATION_DERIVED__

#include "discretization.h"
#include <cfloat>
#include <memory>
#include "operatorAbstract.h"

// forward declaration
class OperatorAbstract;
class IndexFloor;
class OperatorDiagonal;

/** \ingroup Discretizations */

/** @brief Generate at new discretization from a given discretization
 *
 *  A range of transformations allow generation of a new discretization, e.g.:
 *
 *  - grid representation (suitable for plotting) from an original finite element discretization
 *  - representation in terms of the spectral eigenfunctions of some operator
 *  - constrained version of a larger representation
 *
 *  etc.
 *
 *  member operators \b mapFromParent and \b mapToParent map between the discretizations (not necessarily lossless)
 */
class DiscretizationDerived : public Discretization {

protected:
    // ==== data =============================================
    mutable std::shared_ptr<OperatorAbstract> _mapFromParent;
    mutable std::shared_ptr<OperatorAbstract> _mapToParent;

public:

//    /// eigenbasis of Operator select a maximum of MaxN vectors between Emin and Emax
//    DiscretizationDerived(Discretization * Parent, Operator & Ovr, Operator & Op,
//                          double Emin=-DBL_MAX, double Emax=DBL_MAX, int MaxN=INT_MAX, bool excludeEnergyRange = false);

    /// just Pseudo-Orthonormalize given states with respect to projectors and overlap and make new disc from the resulting basis
    DiscretizationDerived(Discretization* Parent, std::vector<Coefficients>& basis);

    virtual ~DiscretizationDerived();

    void setFromParent(std::shared_ptr<OperatorAbstract> Map){_mapFromParent=Map;}
    void setToParent(std::shared_ptr<OperatorAbstract> Map)  {_mapToParent=Map;}

    virtual const OperatorAbstract * mapFromParent() const {return _mapFromParent.get();}
    virtual const OperatorAbstract * mapToParent() const {return _mapToParent.get();} ///< maps back to Parent, not dual or mixed space (i.e. necessary inverses are applied)

protected:
    DiscretizationDerived():_mapFromParent(0),_mapToParent(0){}
public:
};
#endif
