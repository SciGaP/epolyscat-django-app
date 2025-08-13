// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorFloorGP.h"

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

bool OperatorFloorGP::_store = false;
bool OperatorFloorGP::_iterations = false;


/* info
 * nodevals and basis vals are running in SizX which is at _nBeg smaller then the size of _dvrWeights and _dvrPoints
 * therefore while summing over SizX the indices of _dvrWeights and _dvrPoints should be added to _nBeg
 * S^(-1) is multiplied somwhere else in the code for the final axpy
 *
 * Ground Energy 0.869945212029
 */

OperatorFloorGP::OperatorFloorGP(const std::string Name, const std::string Def, const Index *IIndex, const Index *JIndex)
    :OperatorFloor(IIndex->sizeCompute(),JIndex->sizeCompute(),"GrossPitaevskii")
{
}

OperatorFloorGP::OperatorFloorGP(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf)
    : OperatorFloor("GrossPitaevskii")
{
    unpackBasic(Info,Buf);
}

OperatorFloorGP::OperatorFloorGP(std::string Pot, const Index* IIndex, const Index* JIndex, std::complex<double> Multiplier)
    : OperatorFloor("GrossPitaevskii"),
      _nodeVal((dynamic_cast<const BasisDVR*>(IIndex->basis()))->valNodes()),
      _nBeg((dynamic_cast<const BasisDVR*>(IIndex->basis()))->nBeg())
{
    dat=0;
    oNorm=1.;
    // zero norm in absorptive range or when non-local
    if(IIndex->basis()->isAbsorptive()
            or IIndex->basis()->integrable()->lowBound()!=JIndex->basis()->integrable()->lowBound()
            or IIndex->basis()->integrable()->upBound()!=JIndex->basis()->integrable()->upBound()
            )
        oNorm=0.;

    const BasisDVR * bi=dynamic_cast<const BasisDVR*>(IIndex->basis());
    const BasisDVR * bj=dynamic_cast<const BasisDVR*>(JIndex->basis());
    if(bi==0)DEVABORT("not a DVR basis: "+IIndex->basis()->str());
    if(bj==0)DEVABORT("not a DVR basis: "+JIndex->basis()->str());
    if(bi->nodes()!=bj->nodes())DEVABORT("DVR nodes do not match: "+IIndex->basis()->str()+IIndex->basis()->str());
    
    for(std::complex<double> c: bi->valNodes())_basSq.push_back(std::norm(c));
    _GP.assign(bi->size(),0.);
    _GP1.assign(bi->size(),0.);
    
    bi->dvrRule(_dvrPoints,_dvrWeights);
    _bsize=bi->size();
    _idx=IIndex;
}

void OperatorFloorGP::updateNonLin(double time, Coefficients* C){
    Coefficients *Cfloor=C->retrieve(_idx);
    for(int n=0;n<_idx->size();n++){
        _GP1[n]=std::norm(Cfloor->data()[n])*_dvrWeights[n+_nBeg]*_basSq[n]*_basSq[n];
    }
    double norm=0;
    for(const Index* idx=_idx->root()->firstLeaf()->parent();idx!=0;idx=idx->nodeRight()){
        const BasisDVR* bi=dynamic_cast<const BasisDVR*>(idx->basis());
        std::vector<double> dvrPointsI, dvrWeightsI;
        bi->dvrRule(dvrPointsI, dvrWeightsI);
        const Coefficients* Cfloor=C->retrieve(idx);
        for(int j=0; j<Cfloor->size(); j++){
            norm+=std::norm(Cfloor->data()[j])*dvrWeightsI[j+bi->nBeg()]*pow(bi->valNodes()[j],2);
        }
        dvrPointsI.clear();
        dvrWeightsI.clear();
    }
    norm=0;
    for(int j=0; j<Cfloor->size(); j++){
        norm+=std::norm(C->data()[j])*_dvrWeights[j+_nBeg]*pow(_nodeVal[j],2);
    }
}

void OperatorFloorGP::axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX, const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const
{
    if(_iterations){          // Governs whether the iterative method is used
        for(int k=0;k<SizX;k++) Y[k]=Beta*Y[k]+Alfa*_GP1[k]*X[k];
    }
    else {
        for(int k=0;k<SizX;k++){
            Y[k]=Beta*Y[k]+Alfa*_dvrWeights[k+_nBeg]*std::norm(_basSq[k])*std::norm(X[k])*X[k];
        }
    }
}

void OperatorFloorGP::pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const
{
    Buf.insert(Buf.end(),dat->begin(),dat->end());
    packBasic(Info,Buf);
}
