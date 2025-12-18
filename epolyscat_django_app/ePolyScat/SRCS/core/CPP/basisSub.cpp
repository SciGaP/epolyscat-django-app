// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisSub.h"
#include "printOutput.h"
#include "basisIntegrable.h"

BasisSub::BasisSub(const BasisAbstract *Bas, const std::vector<int> &Subset)
    :_bas(Bas)
{
    _name=superBas(Bas)->name();
    if(Subset.size()==0)
        for(unsigned int k=0;k<Bas->size();k++)_subset.push_back(k);
    else{
        _subset=Subset;
    }
    if(_subset.size()==Bas->size()){
        PrintOutput::DEVwarning(Sstr+"identical subset"+Subset,1);
    }
}
unsigned int BasisSub::order() const {
    if(not superBas(this)->integrable())DEVABORT("superBas of BasisSub has no order(): "+superBas(this)->str());
    return superBas(this)->integrable()->order();
}

BasisSub::BasisSub(std::string StrDefinition){
    // StrDefinition format example: "BasisSub: 0,1 of string-definition-of-superbasis"
    if(StrDefinition.find("BasisSub:")!=0)DEVABORT("not a BasisSub, found: "+StrDefinition);
    std::string sub=StrDefinition.substr(StrDefinition.find("BasisSub:")+std::string("BasisSub:").length());
    sub=sub.substr(0,sub.find("of"));
    std::vector<std::string> subStr=tools::splitString(tools::stringInBetween(StrDefinition,"BasisSub:","of",true),',');
    for(std::string s: subStr){
        if(tools::cropString(s)!="")_subset.push_back(tools::string_to_int(s));
    }
    _bas=BasisAbstract::factory(StrDefinition.substr(StrDefinition.find("of")+3));
    _name=_bas->name();
    if(_subset.size()==_bas->size()){
        PrintOutput::DEVwarning(Sstr+"identical subsetA"+_subset,1);
    }
}

std::string BasisSub::strDefinition(const BasisAbstract* Bas, std::vector<int> Subset){
    BasisSub b(Bas,Subset);
    return b.strDefinition();
}

std::string BasisSub::strDefinition() const {
    std::string s("BasisSub: ");
    for(int n: _subset)s+=tools::str(n)+",";
    s=s.substr(0,s.size()-1)+" of "+_bas->strDefinition();
    return s;
}

void BasisSub::valDer(const UseMatrix &X, UseMatrix &Val, UseMatrix &Der, bool ZeroOutside) const{
    UseMatrix val,der;
    if(not _bas->integrable())DEVABORT("cannot evaluate function values of "+_bas->str());
    _bas->integrable()->valDer(X,val,der,ZeroOutside);
    Val=UseMatrix(val.rows(),_subset.size());
    Der=UseMatrix(der.rows(),_subset.size());
    for(size_t k=0;k<_subset.size();k++){
        Val.col(k)=val.col(_subset[k]);
        Der.col(k)=der.col(_subset[k]);
    }
}

bool BasisSub::operator==(const BasisAbstract& Other) const {
    if(this==&Other)return true;
    if(size()!=Other.size())return false;

    const BasisAbstract* sup=superBas(this);
    const BasisAbstract* oth=superBas(&Other);
    if(not sup->operator==(*oth))return false;

    // same super-basis, compare subsets
    const BasisSub* sub=dynamic_cast<const BasisSub*>(&Other);
    if(sub)
        return _subset==sub->_subset;
    else
        for(int s: _subset)if(_subset[s]!=s)return false; // all functions appear in subset in correct ordering
    return true;
}

const BasisAbstract* BasisSub::superBas(const BasisAbstract *Sub){
    const BasisSub* sub=dynamic_cast<const BasisSub*>(Sub);
    if (sub==0)return Sub;
    return superBas(sub->_bas);
}

std::vector<int> BasisSub::subset(const BasisAbstract *Sub){
    const BasisSub* sub=dynamic_cast<const BasisSub*>(Sub);
    std::vector<int> subSet;
    if(sub==0){
        for(size_t k=0;k<Sub->size();k++)subSet.push_back(k);
    }
    else {
        std::vector<int> superSet=subset(superBas(Sub));
        for(int k: sub->_subset)subSet.push_back(superSet[k]);
     }
    return subSet;
}

Eigen::MatrixXcd BasisSub::subMatrix(const Eigen::MatrixXcd &Mat, std::vector<int> ISub, std::vector<int> JSub){
    Eigen::MatrixXcd res(ISub.size(),JSub.size());
    for (size_t j=0;j<JSub.size();j++)
        for (size_t i=0;i<ISub.size();i++)
            res(i,j)=Mat(ISub[i],JSub[j]);
    return res;
}


Eigen::MatrixXcd BasisSub::map(const BasisAbstract *Ibas, const BasisAbstract *Jbas){
    const BasisSub * isub=dynamic_cast<const BasisSub*>(Ibas);
    const BasisSub * jsub=dynamic_cast<const BasisSub*>(Jbas);
    Eigen::MatrixXcd res;
    if((isub==0) == (jsub==0)) return res;
    if(jsub==0){
        if(not(*superBas(isub)==*Jbas))return res;
        res=Eigen::MatrixXcd::Zero(isub->size(),Jbas->size());
        for(size_t k=0;k<isub->size();k++)res(k,isub->_subset[k])=1.;
    }
    else if (isub==0){
        if(not(*superBas(jsub)==*Ibas))return res;
        res=Eigen::MatrixXcd::Zero(Ibas->size(),jsub->size());
        for(size_t k=0;k<jsub->size();k++)res(jsub->_subset[k],k)=1.;
    }
    return res;
}
    
double BasisSub::physical(int Index) const{
    return superBas(this)->physical(subset(this)[Index]);
}
std::string BasisSub::str(int Level) const{
    if(_subset.size()<5)
        return Str("subset{","")+_subset+"} of "+_bas->str(Level);
    else
        return Str("subset{","")+_subset.front()+",...,"+_subset.back()+"} ["+_subset.size()+"] of "+_bas->str(Level);
}
