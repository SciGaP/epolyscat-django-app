// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "volkovGrid.h"

#include "mpiWrapper.h"
#include "discretizationGrid.h"
#include "parameters.h"
#include "coefficientsFunction.h"
#include "wavefunction.h"
#include "operatorIdentity.h"
#include "surfaceFlux.h"
#include "discretizationSurface.h"
#include "operatorMap.h"
#include "basisIntegrable.h"

using namespace std;

VolkovGrid::VolkovGrid(SurfaceFlux *Flux, const Index *SpecIdx):
    OperatorAbstract("Volkov", SpecIdx, SpecIdx){

    string kLevel;
    for(const Index *s=Flux->idx(),*t=SpecIdx;; s=s->child(0), t=t->child(0)){
        if(s->axisName().size()>6 and t->axisName().size()>1){
            if(s->axisName().substr(0,6)=="ValDer" and t->axisName().substr(0,1)=="k"){
                kLevel = t->axisName();
                break;
            }else if(s->axisName().substr(0,6) == "ValDer"){
                if(t->child(0)->axisName().substr(0,1) == "k"){
                    kLevel = t->child(0)->axisName();
                    break;
                }
            }
        }
        if(t->isBottom() or s->isBottom())break;
    }

    // determine axes and quadrature points
    vector<string> axis_names;
    vector<unsigned int> quadPoints;
    for(const Index* ix=Flux->idx();ix!=0;ix=ix->descend()){
        //HACK Not a safe condition, just comparing last character
        if (ix->axisName().find("Phi")!=string::npos or ix->axisName().find("Eta")!=string::npos){
            if(ix->axisName().size()==3 or (ix->axisName().substr(3,ix->axisName().size())==kLevel.substr(kLevel.size()-1,kLevel.size())))
            {
                axis_names.push_back(ix->axisName());
                int maxOrder=0;
                for(const Index *jx=ix;jx!=0;jx=jx->nodeRight()){
                    if(jx->basis()->integrable())maxOrder=std::max(maxOrder,int(jx->basis()->integrable()->order()));
                }
                quadPoints.push_back(maxOrder*2);
            }
        }
    }

    if(axis_names.size()==0){// no axes to be transformed to grid, to enforce equal structure
        gridSpecDisc = new DiscretizationGrid(SpecIdx,{},std::vector<unsigned int>(),{},false);  // create identical discretization
    }
    else{
        gridSpecDisc = new DiscretizationGrid(SpecIdx,axis_names,quadPoints,vector<vector<double>>(),false);
    }

    // Volkov multiplier
    Parameters::update(Flux->FluxBufferBeginTime()); // need time for setup of volkov phase operator
    volkov = new VolkovPhase(gridSpecDisc->idx(),kLevel);

    cGrid = new Coefficients(gridSpecDisc->idx());
    cGrid->treeOrderStorage();

}

void VolkovGrid::update(double Time, const Coefficients* CurrentVec){
    _time = Time;
}

void VolkovGrid::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
        gridSpecDisc->mapFromParent()->apply(1., Vec, 0., *cGrid);
        volkov->multiply(cGrid, cGrid, _time);
        gridSpecDisc->mapToParent()->apply(A, *cGrid, B, Y); // there will be no inverse overlap applied
}

