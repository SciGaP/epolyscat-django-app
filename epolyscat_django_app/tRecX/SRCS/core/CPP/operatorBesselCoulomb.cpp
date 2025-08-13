// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorBesselCoulomb.h"
#include "operatorAbstract.h"
#include "index.h"
#include "basisMat1D.h"
#include "basisBesselCoulomb.h"

OperatorBesselCoulomb::OperatorBesselCoulomb(const std::string& TermOper, const BasisBesselCoulomb* IBas,
                                             const BasisBesselCoulomb* JBas, std::complex<double> Multiplier)
    : OperatorFloor(IBas->size(),JBas->size(),"BesselCoulomb")
{
    BasisMat1D bM1D = BasisMat1D(TermOper,IBas->pure(),JBas->pure());
    std::vector<const Eigen::MatrixXcd*> m = bM1D.mats();
    //...here, possible constraints like banded should be administered...
    std::string hash=TermOper+IBas->hash()+JBas->hash()+"*"+tools::str(Multiplier);
    Eigen::MatrixXcd m0;
    if(m.size()>0){
        if(Multiplier!=1.){
            m0=*m[0]*Multiplier;
            m[0]=&m0;
        }
    }
    else ABORT("m.size()==0 - could not construct operatorBesselCoulomb");
    floor=OperatorFloor::factory(m,hash); // this will recognize and optimally exploit any structure

    // initialize m-indices and b-vectors
    mIdxI = IBas->mIdx(); bVecI = IBas->bVector();
    mIdxJ = JBas->mIdx(); bVecJ = JBas->bVector();

}

void OperatorBesselCoulomb::axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
                                 const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const{
//    if(Alfa!=1.) DEVABORT("Not implemented for Alfa!=1.");
//    if(Beta!=0.) DEVABORT("Not implemented for Beta!=0.");
    if(SizX!=SizY) DEVABORT("SizX!=SizY");

    floor->apply(Alfa,X,SizX,Beta,Y,SizY);
//    // transformation of the j-basis
//    std::vector<std::complex<double> > xTemp;
//    xTemp.assign(SizX,0.);
//    for(unsigned int i=0;i<SizX-1;i++){
//        xTemp[i] = X[i];
//        xTemp[SizX-1] += bVecJ[i]*X[i];
//    }
//    xTemp[SizX-1] += bVecJ[SizX-1]*X[SizX-1];

//    for(unsigned int j=SizX-1;j>mIdxJ;j--){
//        std::swap(xTemp[j],xTemp[j-1]);
//    }

//    // apply the actual floor
//    floor->apply(Alfa,xTemp.data(),SizX,Beta,Y,SizY);

//    // transformation of the i-basis
//    for(unsigned int i=mIdxI;i<SizY-1;i++){
//        std::swap(Y[i],Y[i+1]);
//    }
//    for(unsigned int i=0;i<SizY-1;i++){
//        Y[i] += bVecI[i]*Y[SizY-1];
//    }
//    Y[SizY-1] *= bVecI[SizY-1];
}
