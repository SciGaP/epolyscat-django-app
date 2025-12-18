// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisAssocLeg.h"
#include "basisSetDef.h"


BasisAssocLeg::BasisAssocLeg(std::string Def):BasisIntegrable(-1.,1.)
{
    if(Def.find("AssociatedLegendre:")!=0)DEVABORT("does not define AssociatedLegendre, is "+Def);
    std::vector<std::string> m_size(tools::splitString(Def.substr(std::string("AssociatedLegendre:").length()),','));
    _mAbs=tools::string_to_int(m_size[0]);
    _size=tools::string_to_int(m_size[1]);
}

std::string BasisAssocLeg::strDefinition(const BasisSetDef &Def){
    int mAbs=int(Def.par[0]);
    return std::string("AssociatedLegendre: ")+tools::str(mAbs)+","+tools::str(Def.size());
}

void BasisAssocLeg::valDer(const std::vector<std::complex<double> > &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const{
    OrthogonalNassocLegendre p(_mAbs);
    Val.resize(X.size()*size());
    Der.resize(X.size()*size());
    for(size_t i=0;i<X.size();i++){
        std::vector<double> v,d;
        p.valDer(size(),X[i].real(),v,d);
        for(size_t k=0;k<size();k++){
            Val[i+X.size()*k]=v[k];
            Der[i+X.size()*k]=d[k];
        }
    }
}

void BasisAssocLeg::quadRule(int N, std::vector<double> &QuadX, std::vector<double> &QuadW) const{
    OrthogonalLegendre().quadrature(N,QuadX,QuadW);
}
