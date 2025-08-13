// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "eigenSolverAbstract.h"
#include "eigenSolver.h"
#include "operatorNonLin.h"
#include "operatorFloorGP.h"
#include "operatorMeanEE.h"
#include "operatorDefinition.h"

//#include "operator.h"
#include "basisAbstract.h"
#include "gaunt.h"
#include "coefficientsFloor.h"
#include "index.h"
#include "qtAlglib.h"
//#include "radialmultipole.h"
//#include "basisMat.h"
#include "inverseDvr.h"
#include "basisMat1D.h"
#include "eigenNames.h"
#include "indexNew.h"
#include "basisDvr.h"

class EigenSolverNonLin : public EigenSolverAbstract
{
    void _compute();
    bool *_store;
    bool *_iterations;
    bool *_noInteraction;
    EigenSolver _slv;
    std::string _method;
public:
    EigenSolverNonLin(double Emin, double Emax, int Nmax, bool RightVectors, bool DualVectors, bool ExcludeRange, std::string Method);
    std::vector<std::complex<double> > eigenvalues() {return _slv.eigenvalues();}
    std::vector<Coefficients *> rightVectors() {return _slv.rightVectors();}
    std::vector<Coefficients *> dualVectors() {return _slv.dualVectors();}
};
