// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "discretizationDerived.h"
#include "discretizationGrid.h"
#include "indexDerived.h"
#include "indexGrid.h"
#include "indexNew.h"
#include "operatorFloor.h"
#include "operatorMap.h"
#include "coefficientsFloor.h"
#include "qtEigenDense.h"
#include "printOutput.h"
#include "log.h"

//#include "operator.h" // needed for Old path of discretization
#include "inverse.h"
#include "basisDvr.h"
#include "basisBesselCoulomb.h"
#include "basisSub.h"
#include "basisGridQuad.h"
#include "basisMat1D.h"
#include "operatorWithInverse.h"
#include "timeCritical.h"
#include "mpiWrapper.h"

using namespace std;

DiscretizationGrid::~DiscretizationGrid(){
    for(size_t k=0;k<mapDerivative.size();k++)delete mapDerivative[k];
}

DiscretizationGrid::DiscretizationGrid(const Index *Parent, const std::vector<string> &Axes, std::vector<std::vector<double> > Grid, std::vector<std::vector<double> > Weig, bool Gradient)
    :_parentIndex(Parent){
    _construct(Parent,Axes,Grid,Weig,Gradient);
}


/// convert to grid
///
/// Rules:
///<br> the lowest level with axisName==Axis[k] will be converted to grid
///<br> on each Axis[k], Point[k] equidistant points will be used, if Point[k]==0 default number of coefficients on axis
///<br> Limit[k][l/u=0/1] gives lower an upper bounds for the grid (end-points included), for lb<ub Points[k]>1 needed
///<br> if lb=ub and Points[k] must be > 0, if Points[k]>1, a default interval (axis length) will be used
///<br> Limit.size()==0 and Point.size(), a default quadrature grid will be generated
///<br> Gradient=true will compute gradients wrt to all grid axes
DiscretizationGrid::DiscretizationGrid(const Index *Parent, std::vector<string> Axes, std::vector<unsigned int> Point,
                                       std::vector< std::vector<double> > Boundary, bool Gradient)
    :_parentIndex(Parent)
{
    std::vector<std::vector<double>>grid(Axes.size()),weig(Axes.size());
    IndexGrid::gridWeight(grid,weig,Axes,Point,Boundary); //this should be owned by DiscretizationGrid
    _construct(Parent,Axes,grid,weig,Gradient);
}

void DiscretizationGrid::_construct(const Index *Parent, const std::vector<string> &Axes, std::vector<std::vector<double> > Grid, std::vector<std::vector<double> > Weig, bool Gradient){
    parent=0;
    name=Parent->hierarchy()+"_grid";

    LOG_PUSH("IndexGrid");
    idx()=new IndexGrid(Parent,Axes,Grid,Weig);
    LOG_POP();
    idx()->resetFloor(Parent->firstFloor()->depth()-Parent->depth());

    // purge index
    idx()->purge(Parent->height());
    idx()->sizeCompute();
    // grid maps should always use tensor (not really elegant)
    bool oldTensor=OperatorAbstract::useTensor,oldFuse=OperatorAbstract::fuseOp;
    OperatorAbstract::useTensor=true;
    OperatorAbstract::fuseOp=false;

    _mapFromParent=0;
    _mapToParent=0;
    if(Gradient){
        for(unsigned int k=0;k<Axes.size();k++)
            mapDerivative.push_back(new OperatorMap(idx(),Parent,Axes[k]));
    }

    // switch back to previous settings
    OperatorAbstract::useTensor=oldTensor;
    OperatorAbstract::fuseOp=oldFuse;
    if(Gradient){
        for(unsigned int k=0;k<Axes.size();k++)
            mapDerivative.push_back(new OperatorMap(idx(),Parent,Axes[k]));
    }

}

const OperatorAbstract *DiscretizationGrid::mapFromParent() const {
    if(_mapFromParent==0){
        timeCritical::suspend();
        _mapFromParent.reset(new OperatorMap(idx(),_parentIndex));
        timeCritical::resume();
    }
    return _mapFromParent.get();
}

bool isMapToParent(const Index* Grid, const Index* Parent){
    if(not (Parent->basis()==Grid->basis())){
        if(Parent->basis()->integrable() and
                not BasisMat1D("1",Parent->basis()->integrable(),Parent->basis()->integrable()).mat().isIdentity(1.e-12))return false;
        if(Grid->isLeaf()!=Parent->isLeaf())return false;
        if(not Grid->isLeaf() and not isMapToParent(Grid->child(0),Parent->child(0)))return false;
    }
    else
        for(size_t k=0;k<Grid->childSize();k++)
            if(not isMapToParent(Grid->child(k),Parent->child(k)))return false;

    return true;
}

const OperatorAbstract* DiscretizationGrid::mapToParent() const {
    if(_mapToParent==0){
        timeCritical::suspend();
        if(isMapToParent(idx(),_parentIndex))
            _mapToParent.reset(new OperatorMap(_parentIndex,idx()));
        else {
            _mapToDual.reset(new  OperatorMap(_parentIndex,idx()));
            _mapToParent.reset(new OperatorWithInverse(_mapToDual));
        }
        timeCritical::resume();
    }
    return _mapToParent.get();
}











