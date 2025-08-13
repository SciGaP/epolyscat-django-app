// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisDvrSqrt.h"

#include "printOutput.h"
#include "useMatrix.h"
#include "basisSetDef.h"
#include "tRecXchecks.h"

using namespace std;

BasisDVRSqrt::BasisDVRSqrt(const BasisSetDef &Def):BasisDVR(Def) {
}

void BasisDVRSqrt::valDer(const std::vector<std::complex<double> >&X,
                      std::vector<std::complex<double> >&Val,std::vector<std::complex<double> >&Der, bool ZeroOutside) const{
    // Factor(x)*poly(x)*sqrt(measure(x))

    // starting values sqrt(v(q(x)))*Factor(x), d/dx ( sqrt(v(q(x))*Factor(x) )
    std::vector<double> val0,der0;
    for(size_t i=0;i<X.size();i++){
        double vFactor = sqrt(X[i].real());
        double dFactor = 0.5/sqrt(X[i].real());
        
        double vWeight = _opol->weight(   qFromX(X[i].real()));
        double dWeight = _opol->derWeight(qFromX(X[i].real())) / _scale;

        val0.push_back(sqrt(vWeight)             * vFactor);
        der0.push_back(0.5/sqrt(vWeight)*dWeight * vFactor
                    +  sqrt(vWeight)             * dFactor);
    }

    Val.clear();
    Der.clear();
    // column[n]...function numbers, row[i]...x-values
    for(int n=_nBeg;n<_nBeg+_size;n++){
        // tabulate divisions
        std::vector<double>qDiffN(_dvrX.size(),0.);
        for(size_t k=0;k<_dvrX.size();k++)
            if(k!=size_t(n))qDiffN[k]=1./(_dvrX[n]-_dvrX[k]);

        // get norm of n'th function; on margin, force value=1, else ||b_n||=1
        double vWeight = _opol->weight(   qFromX(_dvrX[n]));
        double vFactor = sqrt(_dvrX[n]);
        
        double normN=1./(vWeight*vFactor*vFactor);
        if(not (   abs(_lowBound-_dvrX[n])<abs(_scale)*1e-10
                   or abs(_upBound-_dvrX[n])<abs(_scale)*1e-10))normN/=_dvrW[n];
        normN=sqrt(normN);

        for(size_t i=0;i<X.size();i++){
            if(ZeroOutside and (X[i].real()<_lowBound or X[i].real()>_upBound)){
                Val.push_back(0.);
                Der.push_back(0.);
            }
            else {
                Val.push_back(val0[i]*normN);
                Der.push_back(der0[i]*normN);
                for(size_t k=0;k<_dvrX.size();k++){
                    if(k==size_t(n))continue;
                    Der.back()=(Der.back()*(X[i].real()-_dvrX[k])+Val.back())*qDiffN[k];
                    Val.back()=(Val.back()*(X[i].real()-_dvrX[k])           )*qDiffN[k];
                }
            }
        }
    }
}

std::string BasisDVRSqrt::str(int Level) const{
    return Str("DVRSqrt(","")+_opol->name()+") ["+_lowBound+","+_upBound+"] ("+_shiftX+","+_scale+") "+_size+"["+_dvrX.size()+"]";
}

bool BasisDVRSqrt::operator==(const BasisAbstract & Other) const {
    if(this==&Other)return true;
    const BasisDVRSqrt* other;
    if(0==(other=dynamic_cast<const BasisDVRSqrt*>(&Other)))return false;
    if(_dvrX!=other->_dvrX)return false;
    if(_dvrW!=other->_dvrW)return false;
    return _opol->name()==other->_opol->name();
}
