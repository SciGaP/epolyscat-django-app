// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "densityOfStates.h"

#include "basisAbstract.h"
#include "basisIntegrable.h"
#include "basisAbstract.h"

#include "operatorAbstract.h"
#include "eigenSolver.h"
#include "readInput.h"
#include "stringTools.h"
#include "asciiFile.h"
#include "printOutput.h"

DensityOfStates::DensityOfStates(ReadInput &Inp):nPoints(0){
    if(not Inp.found("DensityOfStates"))return;
    Inp.read("DensityOfStates","eMin",eMin,"0","return DOS up to eMax");
    Inp.read("DensityOfStates","eMax",eMax,tools::str(DBL_MAX),"return DOS up to eMax");
    Inp.read("DensityOfStates","nPoints",nPoints,"100","number of energy points in (0,eMax]");

    outDir=Inp.outputTopDir();
}

std::vector<double> DensityOfStates::compute(const OperatorAbstract *Op){
    EigenSolver slv(-DBL_MAX,DBL_MAX,true,false);
    slv.compute(Op);
    slv.sort("SmallReal");
    std::vector<std::complex<double> > eval;
    eval=slv.eigenvalues();

    std::vector<std::string> com;
    com.push_back("Eigenvalues of "+Op->name);
    std::vector<std::vector<double> >cols(2);
    for(size_t k=0;k<eval.size();k++){
        cols[0].push_back(std::real(eval[k]));
        cols[1].push_back(std::imag(eval[k]));
    }
    AsciiFile f(outDir+"eig");
    f.writeComments(com);
    f.writeCols(cols);

    box=1.; //HACK for box size
    for(const Index* idx=Op->iIndex;not idx->isLeaf();idx=idx->descend()){
        DEVABORT("HACK needs fixing");
//        if(idx->basisSet()!=0 and idx->basisSet()->eta()!=1.)box*=idx->basisSet()->upBound();
    }

    // up to energy where 90 % of eigenvalues are below
    if(eMax>DBL_MAX/2)eMax=real(eval[int(eval.size()*0.9)]);

        std::vector<double> dens;
        double deltaE=eMax/nPoints;
        for(int k=0;k<nPoints;k++){
            dens.push_back(0.);
            for(size_t l=0;l<eval.size();l++)
                dens.back()+=std::norm(1./(box*(eval[l]-eMin-(k+1)*deltaE)));
        }
    return dens;
}

void DensityOfStates::output(const OperatorAbstract *Op){
    if(nPoints==0)return; // do not compute DOS

    AsciiFile f(outDir+"dos");
    std::vector<std::string> com;
    com.push_back("Density of states for "+Op->name);
    f.writeComments(com);

    std::vector<std::vector<double> >cols(2);
    for(int k=0;k<nPoints;k++)cols[0].push_back(eMin+(k+1)*((eMax-eMin)/nPoints));
    cols[1]=compute(Op);
    f.writeCols(cols);

//    std::vector<std::vector<double> >eval(2);


    PrintOutput::title("density-of-states on file "+f.name());
    exit(0);
}
