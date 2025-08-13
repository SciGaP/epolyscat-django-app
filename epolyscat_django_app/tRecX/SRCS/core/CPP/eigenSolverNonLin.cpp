// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "eigenSolverNonLin.h"
#include "basisGrid.h"
#include "chrono"

/* INFO
 * 1) for norm use psi.scalarProduct(psi), b.c. norm() is for propagation and innerProduct is just coefficients-wise
 * 2) to normalise scale with 1./sqrt(psi.scalarProduct(psi))
 * 3) qe-q is 0.625
 * 3) qe-qcos is 0.66668512
 *
 */

EigenSolverNonLin::EigenSolverNonLin (double Emin, double Emax, int Nmax, bool RightVectors, bool DualVectors, bool ExcludeRange, std::string Method)
    : _slv(Emin, Emax, Nmax, RightVectors, DualVectors, ExcludeRange, Method)
{
}

void EigenSolverNonLin::_compute()
{
    //auto start = std::chrono::steady_clock::now();

    const OperatorTree* OpTree=dynamic_cast<const OperatorTree*>(_op);
    const OperatorTree* OvrTree=dynamic_cast<const OperatorTree*>(_ovr);

    const OperatorTree OpZ("Z","<<Z>>",_op->iIndex,_op->jIndex);

    Coefficients psi(_op->jIndex);
    psi.setToZero();



    std::complex<double> E=0;
    std::complex<double> Eprev;
    std::string chooseOperator="EE";
    if (chooseOperator=="GP"){
        DEVABORT("reimplement update for propagation as in MeanEE");
        _iterations=&OperatorFloorGP::_iterations;
    }
    if (chooseOperator=="EE"){
        _noInteraction=&OperatorMeanEE::_noInteraction;
    }

    //-------------------------------------------------------------------------------------------------


    //------------------------------------------------------------------------------------------------------------------------------

    std::complex<double> matrixElement=0.;
    int i=0;
    OperatorMeanEE::_iterations=true;
    do{
        i++;
        _slv.clear();
        if(false){

        }
        else{
            Eprev=E.real();
            if(chooseOperator=="EE"){
                *_noInteraction=true;
                matrixElement=OpTree->matrixElement(psi,psi);
                *_noInteraction=false;
                OperatorMeanEE::ME=matrixElement;
            }
            const_cast<OperatorTree*>(OpTree)->updateNonLin(0., &psi);
            _slv.compute(OpTree,OvrTree);
            _slv.normalize();
            //OperatorMeanEE::ME=matrixElement=0;
        }

        _slv.select("SmallReal[1]");
        psi=*_slv.rightVectors()[0];

        _slv.eigenvalues()[0];
        E=_slv.eigenvalues()[0];

        if(i==200)DEVABORT("Energy does not converge, i=200");
    }while(std::abs(E.real()-Eprev)>1.e-2);
    OperatorMeanEE::_iterations=false;
    OperatorMeanEE::ME=matrixElement=0;
    _eigenvalues=_slv.eigenvalues();
    _rightVectors=_slv.rightVectors();
    _leftVectors=_slv.rightVectors();
    _dualVectors=_slv.dualVectors();

    select("SmallReal[1]");
    _eigenvalues[0]=_slv.eigenvalues()[0]+matrixElement;


    //psi.setToZero();
    //const_cast<OperatorTree*>(OpTree)->updateNonLin(0., &psi);
}
