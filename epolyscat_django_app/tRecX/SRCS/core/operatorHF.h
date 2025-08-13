// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORHF_H
#define OPERATORHF_H

#include <memory>
#include "operatorAbstract.h"
#include "operatorNonLin.h"
#include "operatorFloorNonLin.h"
#include "operatorHartree.h"
#include "operatorDefinition.h"

class OperatorFloorHF : public OperatorFloorNonLin
{
    // for now, "has a" OperatorHartree
    static std::shared_ptr<OperatorHartree> _hfOp;
public:
    /// single Hartree-Fock floor (needs postProcess)
    OperatorFloorHF(std::string Pot, const Index *IIndex, const Index *JIndex,std::complex<double> Multiplier);
    void updateNonLin(double Time, const Coefficients* C);

    void axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
              const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const{DEVABORT("not implemented");}
    void pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const{DEVABORT("not implemented");}

    static void postProcess(OperatorTree * Op);

    /// hartree matrix element for product function (slow, for debug)
    static std::complex<double> matrixElement(const Coefficients* C);
};

class OperatorHF:public OperatorNonLin {
    std::shared_ptr<OperatorTree> _hf;
    std::shared_ptr<Coefficients> _cChan;
    std::shared_ptr<Index> _iChan;
    OperatorDefinition _def;
    std::shared_ptr<Coefficients> _applyV,_applyY;
public:
    OperatorHF(std::string Pot, const Index *IIndex, const Index *JIndex);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    void update(double Time, const Coefficients *CurrentVec);
};

#endif // OPERATORHF_H
