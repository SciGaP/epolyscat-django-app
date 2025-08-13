// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "multipoleOrbital.h"
#include "coefficients.h"
#include "basisOrbitalNumerical.h"
#include "basisExpIm.h"
#include "index.h"
std::map<std::string,std::vector<std::map<int,std::vector<int>>>> MultipoleOrbital::_listLangle;
std::map<std::string,Eigen::MatrixXcd> MultipoleOrbital::_listRho;
std::map<std::string,std::vector<std::map<int,std::vector<Coefficients*>>>> MultipoleOrbital::_listCoef;

int MultipoleOrbital::channelIndex(const std::string&Chan,const Index* Idx){
    if(Idx->parent()==0)ABORT("channel level "+Chan+"not foun in"+Idx->hierarchy());
    if(Idx->parent()->axisName()==Chan)return Idx->nSibling();
    return channelIndex(Chan,Idx->parent());
}

std::complex<double> MultipoleOrbital::density(std::string Name,const Index* IIndex, const Index* JIndex) const{
    return _rho->operator()(channelIndex(Name,IIndex),channelIndex(Name,JIndex));
}

int MultipoleOrbital::mQuantumNumber(const Index* Idx){
    if(Idx->parent()==0)ABORT("no Phi-axis found in "+Idx->hierarchy());
    if(Idx->parent()->axisName()=="Phi"){
        const BasisExpIm* b=dynamic_cast<const BasisExpIm*>(Idx->parent()->basis());
        if(not b)ABORT("mQuantumNumber only for BasisExpIm, is: "+Idx->parent()->basis()->str());
        return b->physical(Idx->nSibling());
    }
    return mQuantumNumber(Idx->parent());
}

int MultipoleOrbital::lQuantumNumber(const Index* Idx){
    if(Idx->parent()==0)ABORT("no Eta-axis found in "+Idx->hierarchy());
    if(Idx->parent()->axisName()=="Eta"){
        return Idx->parent()->basis()->physical(Idx->nSibling());
    }
    return lQuantumNumber(Idx->parent());
}

void MultipoleOrbital::add(std::string Name, const BasisOrbitalNumerical &Orbs, const std::vector<std::complex<double> > Rho){
    if(_listCoef.count(Name))DEVABORT("multiple use of name for MultipoleOrbital "+Name);
    _listCoef[Name].resize(Orbs.size());
    _listLangle[Name].resize(Orbs.size());
    if(Rho.size()==0)_listRho[Name]=Eigen::MatrixXcd::Identity(Orbs.size(),Orbs.size());
    else {
        if(Rho.size()!=std::pow(Orbs.size(),2))
            ABORT(Str("Rho-matrix size does not match orbitals: ")+Rho.size()+"!="+Orbs.size()+"x"+Orbs.size());
        _listRho[Name]=Eigen::Map<Eigen::MatrixXcd>(const_cast<std::complex<double>*>(Rho.data()),Orbs.size(),Orbs.size());
    }
    for(size_t i=0;i<Orbs.size();i++){
        Coefficients* c=const_cast<Coefficients*>(Orbs.orbital(i)->firstLeaf());
        for(;c!=0;c=c->nextLeaf()){
            int m=mQuantumNumber(c->idx());
        _listCoef[Name][i][m].push_back(c);
        _listLangle[Name][i][m].push_back(lQuantumNumber(c->idx()));
        }
    }
}

MultipoleOrbital::MultipoleOrbital(std::string Name, const Index *Idx)
{
    if(not _listCoef.count(Name))
        ABORT(Name+" not in MultipoleOrbital, available: "+tools::listMapKeys(_listCoef));

    _lCoef=&_listCoef[Name][channelIndex(Name,Idx)][mQuantumNumber(Idx)];
    _lAngle=&_listLangle[Name][channelIndex(Name,Idx)][mQuantumNumber(Idx)];
    _rho=&_listRho[Name];
}


//const Eigen::MatrixXcd MultipoleOrbital::vals(int K) const {
//    return Eigen::Map<Eigen::MatrixXcd>(_lCoef->data()[K]->data(),_lCoef->data()[K]->size(),1);
//}
