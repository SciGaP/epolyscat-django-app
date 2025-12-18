// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include <fstream>

#include "readInput.h"
#include "printOutput.h"
#include "units.h"
#include "mpiWrapper.h"

#include "operatorTree.h"
#include "operatorRALL.h"
#include "index.h"
#include "eigenTools.h"

using namespace std;

int main(int argc, char* argv[]) {
    MPIwrapper::Init(argc,argv);

    //============================================================================================
    // input and part of setup
    //==============================================================================================

    //     ReadInput::openMain("",argc,argv);
    //    Units::setDefault("au");        // general default units
    //    ReadInput::main.setUnits("au"); // convert input to these units

    std::string name("CACHE_TRECX/small_hamiltonian_Rn.Phi.Eta.Rn_Rn.Phi.Eta.Rn");
    name="CACHE_TRECX/rall_hamiltonian_Rn.Phi.Eta.Rn_Rn.Phi.Eta.Rn";
    std::ifstream opFil(name,std::ios_base::binary);
    OperatorTree op(opFil);
    opFil.close();
    opFil.open(name,std::ios_base::binary);
    OperatorTree opc(opFil,op.name,op.iIndex,op.jIndex);
    OperatorRALL *rall=new OperatorRALL(op,{0},1,50);

    if(rall->iIndex->size()<1000){
        Eigen::MatrixXcd mr=rall->matrix();
        Eigen::MatrixXcd m0=opc.matrix();
        if((mr-m0).isZero(1.e-12))
            Sout+"OK"+Sendl;
        else
            Sout+"fail"+EigenTools::str(m0,1)+Sendl;
    }

    Coefficients cx(op.jIndex),cy(op.iIndex),ccy(op.iIndex);
    cx.setToRandom();
    rall->apply(1,cx,0,cy);
    opc.apply(1,cx,1,ccy);
    cy-=ccy;
    cy=cy.cwiseDivide(ccy);
    if(cy.isZero(1.e-6))
        Sout+"fine"+Sendl;
    else
        Sout+cy.str(0)+Sendl;

    Sout+"all done"+Sendl;
}
