// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONGRID_H
#define DISCRETIZATIONGRID_H

#include "discretizationDerived.h"
#include "indexDerived.h"
#include "operatorAbstract.h"
#include "tree.h"
#include <memory>

class OperatorFloor;
class OperatorIdentity;

/** \ingroup Discretizations */
/// Discretization by converting function on selected axes to grids
class DiscretizationGrid:public DiscretizationDerived
{
    const Index* _parentIndex;
    mutable std::shared_ptr<const OperatorAbstract> _mapToDual; // mutable as it will be created upon first use
    void _construct(const Index *Parent, const std::vector<std::string> &Axes, std::vector<std::vector<double>> Grid, std::vector<std::vector<double>> Weig, bool Gradient);
public:

    virtual ~DiscretizationGrid();
    DiscretizationGrid(){}

    /// standard constructor
    DiscretizationGrid(const Index *Parent, const std::vector<std::string> &Axes, std::vector<std::vector<double>> Grid, std::vector<std::vector<double>> Weig, bool Gradient);

    /// alternate constructor, use Index for parent (legacy)
    DiscretizationGrid(const Index *Parent /** original discretizaiton */,
                       std::vector<std::string> Axes /** names of axes to be transformed */,
                       std::vector<unsigned int> Point=std::vector<unsigned int>() /** number of points, minimal quadrature grid if not specified */,
                       std::vector<std::vector<double> > Limit=std::vector<std::vector<double> >() /** lower and upper grid boundaries */,
                       bool Gradient=false /** also compute transformations to gradients wrt all grid axes */);

    const OperatorAbstract *mapFromParent() const;
    const OperatorAbstract *mapToParent() const; ///< maps back to one of To = Parent, Dual, Mix

    std::vector<const OperatorAbstract*> mapDerivative;
};

#endif // DISCRETIZATIONGRID_H
