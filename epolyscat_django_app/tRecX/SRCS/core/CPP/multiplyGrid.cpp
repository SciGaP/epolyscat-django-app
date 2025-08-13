// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "multiplyGrid.h"
#include <vector>

#include "discretizationDerived.h"
//#include "operator.h"
#include "coefficients.h"
#include "parameters.h"
#include "coefficientsFloor.h"
#include "coefficientsFunction.h"
#include "index.h"
#include "discretizationGrid.h"
#include "printOutput.h"

using namespace std;

MultiplyGrid::MultiplyGrid(const Discretization *D /** original, non-grid discretization */,
                           std::string Func /** string specifying the function (see "function factory" in code) */,
                           std::string Grid /** axis to transform: format as axisName1.axisName2.axisName3 */):
    OperatorAbstract(Func,D->idx(),D->idx())
{
    // get parameters for grid transformation
    vector<string>level=tools::splitString(Grid,'.');
    vector<string>axes=tools::splitString(D->name,'.');
    vector<int> transPar(2*level.size(),-1);
    for(unsigned int l=0;l<level.size();l++){
        const Index* I=D->idx();
        for(unsigned int k=0;k<axes.size();k++){
            if(axes[k]==level[l]){
                transPar[2*l]=k;
                transPar[2*l+1]=0;
            }
            I=I->child(0);
        }
        if(transPar[2*l]==-1)ABORT("could not find axis "+level[l]+" in discretization "+D->name);
    }

    //    dGrid=new DiscretizationDerived(D,DiscretizationDerived::grid,transPar);
    PrintOutput::warning("grid transformation in MultiplyGrid needs to be checked after transition to DiscretizationGrid");
    dGrid=new DiscretizationGrid(D->idx(),level);
    dGrid->mapFromParent()->lhsVector();

    // function factory
    if(     Func=="Identity")       function=new Identity();  // for testing
    else ABORT("function not defined: "+Func);
}

void MultiplyGrid::apply(std::complex<double> Alfa, const Coefficients &X, std::complex<double> Beta, Coefficients &Y) const
{
    if(dGrid==0)ABORT("not implemented");
    dGrid->mapFromParent()->apply(1.,X,0.,*dGrid->mapFromParent()->tempLHS());
    function->multiply(dGrid->mapFromParent()->tempLHS(),dGrid->mapFromParent()->tempLHS(),_time);
    dGrid->mapToParent()->apply(Alfa,*dGrid->mapFromParent()->tempLHS(),Beta,Y);
}

