// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorVectors.h"

#include "index.h"
#include "coefficients.h"
#include "timer.h"
OperatorVectors::OperatorVectors(std::string Name, const Index *IIndex, const Index *JIndex)
    :OperatorAbstract(Name,IIndex,JIndex){vecs.resize(JIndex->sizeCompute(),Coefficients(IIndex,0.));}

void OperatorVectors::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    Y.scale(B);
    for(int k=0;k<vecs.size();k++)
        Y.axpy(A*const_cast<Coefficients*>(&Vec)->data()[k],&vecs[k]);
}

void OperatorVectors::insertColumn(unsigned int J, const Coefficients  &C){
    if(C.idx()!=iIndex)ABORT("cannot insert, index does not match");
    vecs[J]=C;
}
