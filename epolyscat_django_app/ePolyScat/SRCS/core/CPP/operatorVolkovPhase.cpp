// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorVolkovPhase.h"
#include "discretizationDerived.h"
#include "coefficientsFunction.h"
#include "parameters.h"

OperatorVolkovPhase::OperatorVolkovPhase(const DiscretizationDerived* GridSpecDisc, const std::string KLevel):OperatorAbstract("OperatorVolkovPhase",GridSpecDisc->idx(),GridSpecDisc->idx()){
    volkov = new VolkovPhase(GridSpecDisc,KLevel);
}

void OperatorVolkovPhase::apply(const std::complex<double> A, const Coefficients& X, const std::complex<double> B, Coefficients& Y) const {
    volkov->multiply(&X,&Y,_time);
}

void OperatorVolkovPhase::update(double Time, const Coefficients* CurrentVec){
    _time=Time;
    Parameters::update(Time);
}
