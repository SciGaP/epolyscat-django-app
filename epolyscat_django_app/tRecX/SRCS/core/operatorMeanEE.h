// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorFloorNonLin.h"
#include "basicDisc.h"
#include "readInput.h"
#include "useMatrix.h"
#include "operatorZD.h"

//#include "operator.h"
#include "basisAbstract.h"
#include "gaunts.h"
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
#include "testIntegration.h"
#include "multipolePotential.h"
#include "basisAssocLeg.h"

class OperatorMeanEE : public OperatorFloorNonLin
{
    void computeMatMultiClass(const Index* idx); //Using more general MultipolePotential class
    void computeRadialMatrix(const Index *IIndex, const Index *JIndex);

    std::vector<std::vector<std::vector<double> > > _R;            // Radial matrix Rij=\int\int bibj/|r-r|
    std::vector<std::complex<double>> _W;           //diagonal matrix for axpy-application
    double _matrixElement;                         // mean energy <psi|H0|psi> for the first electron without interaction

    const BasisDVR* _b; //current floor basis
    std::vector<std::vector<std::vector<double>>> _GauntMatrix;

    std::complex<double> getC(const Coefficients* C, int i, int j); // get an element C_lj from coefficient
    void getMeanField(const Coefficients *C); // compute _W

    int _etaJ; // current l of the floor wrt J
    int _etaI; // current l of the floor wrt I
    int _lsize; // order of Eta

    std::vector<std::complex<double> >_ovr; //needed for speed up
public:
    OperatorMeanEE(const std::string Name, const std::string Def, const Index *IIndex, const Index *JIndex);
    OperatorMeanEE(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf);
    OperatorMeanEE(std::string Pot, const Index* IIndex, const Index* JIndex, std::complex<double> Multiplier);

    void updateNonLin(double Time, const Coefficients *C);
    void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
              const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const;

    void pack(std::vector<int> & Info,std::vector<std::complex<double> > &Buf) const;

    //Switchers
    static bool disableUpdate;
    static bool _iterations;
    static bool _noInteraction;
    static std::complex<double> ME;

    void test();
};

