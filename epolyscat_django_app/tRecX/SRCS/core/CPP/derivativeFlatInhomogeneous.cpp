// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "derivativeFlatInhomogeneous.h"
#include "tsurffSource.h"
//#include "operator.h"
#include "timer.h"
#include "printOutput.h"

#include "inverse.h" //DEBUG only
#include "debugInfo.h"

using namespace std;

DerivativeFlatInhomogeneous::~DerivativeFlatInhomogeneous(){
    for(auto sf: tS)delete sf;
}

DerivativeFlatInhomogeneous::DerivativeFlatInhomogeneous(const OperatorTree *Op, double ApplyThreshold, const DiscretizationSpectral *ProjectionDisc, std::vector<TsurffSource*> Source):
    tS(Source),_op(Op)
{
    applyEpsilon=ApplyThreshold;
    if(Op!=0)_construct(Op,ProjectionDisc,0);
    else if(tS.size()!=0)PrintOutput::message("Only integration of source - no Propagation operator defined");
    else DEVABORT("neither source nor operator defined - cannot do anything");
}

DerivativeFlatInhomogeneous::DerivativeFlatInhomogeneous(const OperatorTree *Op, double ApplyThreshold,
                                                         std::shared_ptr<ProjectSubspace> Projection, std::vector<TsurffSource*> Source):
    tS(Source),_op(Op)
{
    applyEpsilon=ApplyThreshold;
    if(Op!=0)_construct(Op,0,Projection);
    else if(tS.size()!=0)PrintOutput::message("Only integration of source - no Propagation operator defined");
    else DEVABORT("neither source nor operator defined - cannot do anything");
}

void DerivativeFlatInhomogeneous::update(double Time, const Coefficients* CurrentVec)
{
    DerivativeFlat::update(Time,CurrentVec);
    // Update Source
    for(unsigned int i=0;i<tS.size();i++)
        if(tS[i]->UpdateSource(Time)==0)
            ABORT("Could not update Source");
}

void DerivativeFlatInhomogeneous::apply(std::complex<double> A, const Coefficients &X, std::complex<double> B, Coefficients &Y) const
{
    DEVABORT("who's calling?");
}

void DerivativeFlatInhomogeneous::apply(std::complex<double> A, CoefficientsLocal *localX, std::complex<double> B, CoefficientsLocal &Y) const
{
    if(not isEmpty()){
        DerivativeFlat::applyA(A,localX);
        for(unsigned int i=0;i<tS.size();i++)localXY->axpy(A,tS[i]->CurrentSource());
        DerivativeFlat::applyB(B,Y);

    }
    else {
        Y.scale(B);
        for(unsigned int i=0;i<tS.size();i++)Y.axpy(A,tS[i]->CurrentSource());
    }
}

const Index* DerivativeFlatInhomogeneous::idx() const {
    return  DerivativeFlat::idx() ? DerivativeFlat::idx() : tS[0]->idx();
}
