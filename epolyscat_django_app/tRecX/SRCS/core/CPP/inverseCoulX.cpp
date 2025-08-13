// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "inverseCoulX.h"
#include "index.h"
#include "operatorInverseCoulX.h"

InverseCoulX::InverseCoulX(const Index *Idx, unsigned int SubD, unsigned int SuperD, bool BandOvr):Inverse("Inv("+Idx->hierarchy()+")",Idx,Idx){
    opInvCoulX = new OperatorInverseCoulX(Idx,SubD,SuperD,BandOvr);
}

void InverseCoulX::apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    opInvCoulX->apply(A,Vec,B,Y);
}
