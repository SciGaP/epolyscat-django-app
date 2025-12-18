// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorNonlinPsi.h"

#include "index.h"
#include "basisAbstract.h"
#include "operatorMap.h"

OperatorNonlinPsi::OperatorNonlinPsi(std::string Def, const Index *IIndex, const Index *JIndex)
    :OperatorFloor(IIndex->sizeStored(),JIndex->sizeStored(),"NonlinPsi")
{
    if(IIndex!=JIndex)DEVABORT("OperatorNonlinPsi is strictly local");
    if(IIndex->basis()->grid()==0 or JIndex->basis()->grid()==0)DEVABORT("OperatorNonlinPsi only for grid (at present)");
}
void OperatorNonlinPsi::axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
                             const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const{
    for(int k=0;k<SizX;k++)Y[k]=Beta*Y[k]+Alfa*X[k]*std::norm(X[k]);
}
