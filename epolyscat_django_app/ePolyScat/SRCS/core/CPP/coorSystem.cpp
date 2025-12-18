// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "coorSystem.h"

#include "tools.h"
#include "str.h"
#include "coordinate.h"
#include "coorPolar3D.h"
#include "coorCylinder3D.h"
#include "coorRhoZ.h"
#include "coorEtaR.h"
#include "coorEtaR_PEtaPXi.h"

using namespace std;

static vector<string> sysAvail;
CoorSystem* CoorSystem::trySystem(CoorSystem* & S){
    if(find(sysAvail.begin(),sysAvail.end(),S->_standardName)==sysAvail.end())sysAvail.push_back(S->_standardName);
    if(S->posInStandard().size()!=0)return S;
    delete S;
    S=0;
    return 0;
}

CoorSystem* CoorSystem::factory(std::string System){
    CoorSystem * s;
    if(trySystem(s=new CoorPolar3D(System)))return s;
    if(trySystem(s=new CoorCylinder3D(System)))return s;
    if(trySystem(s=new CoorEtaR_PEtaPXi(System)))return s;
    if(trySystem(s=new CoorEtaR(System)))return s;
    if(trySystem(s=new CoorRhoZ(System)))return s;
    ABORT(Str(("Coordinate system not defined: "+System),"\n")+"must be (numbered) permutation of any of: "+sysAvail);
    return 0;
}

CoorSystem* CoorSystem::factory(std::string OutSys, std::string InSys){
    // primitive version, can be improved
    CoorSystem * s=factory(OutSys);
    if(not s or s->refSystem()!=InSys)
        ABORT(Str(("Coordinates "+OutSys+" not reacheable from "+InSys),"\n")+"must be (numbered) permutation of any of: "+sysAvail);
    return s;
}

CoorSystem::CoorSystem(std::string StandardName, std::string RefCoor, string Name):_standardName(StandardName),_name(Name),_ref(RefCoor)
{
    if(_name=="")_name=_standardName;
}

static std::map<std::string,std::string> aka={
    {"Rn","R"}
};
vector<string> CoorSystem::standardCoors(std::string Name) const{

    vector<string> stdC=tools::splitString(Name,'.');
    string num;
    for(auto& c: stdC){
        // strip number
        size_t p;
        if((p=c.find_first_of("0123456789"))!=string::npos){
            if(p!=c.length()-1)ABORT("coordinates can have at most a single digit at end, found: "+c);
            if(num!="" and num!=c.substr(p))ABORT("all coordinates in system must have same number, is: "+Name);
            num=c.substr(p);
            c=c.substr(0,p);
        }
        // substitute alias
        c=Coordinate::kind(c);
        if(aka.count(c))c=aka[c];
    }
    return stdC;
}

std::vector<int> CoorSystem::posInStandard() const {
    vector<string> namS=standardCoors(_name);
    vector<string> stdS=tools::splitString(_standardName,'.');
    if(namS.size()!=stdS.size())return vector<int>();

    vector<int> pos;
    for(size_t k=0,p;k<stdS.size();k++){
        pos.push_back((p=std::find(stdS.begin(),stdS.end(),namS[k])-stdS.begin()));
        if(p==stdS.size())return vector<int>();
    }
    return pos;
}

vector<double> CoorSystem::fromRef(const std::vector<double> &Ref) const{
    vector<int> pos=posInStandard();

    vector<double> stdCoor=_fromRef(Ref);
    vector<double> coor(pos.size());
    for (size_t k=0;k<pos.size();k++)coor[k]=stdCoor[pos[k]];
    return coor;
}

vector<double> CoorSystem::toRef(const std::vector<double> &Coor) const{
    vector<int> pos=posInStandard();

    vector<double> stdCoor(pos.size());
    for (size_t k=0;k<pos.size();k++)stdCoor[pos[k]]=Coor[k];
    return _toRef(stdCoor);
}

Eigen::MatrixXd CoorSystem::jacRefdCoor(const std::vector<double> &Coor) const{
    vector<int> pos=posInStandard();

    vector<double> stdCoor(pos.size());
    for (size_t k=0;k<pos.size();k++)stdCoor[pos[k]]=Coor[k];

    vector<double> stdJac=_jacRefdCoor(stdCoor);

    // standard indexing in jacobian: J[ij]=(df/dq)_ij = df_i/dq_j
    // sort columns from standard to actual sorting of q_j
    vector<double> jac(stdJac.size());
    for(size_t j=0,ij=0;j<pos.size();j++)
        for (size_t i=0;i<pos.size();i++,ij++){
            jac[ij]=stdJac[pos.size()*pos[j]+i];
        }
    return Eigen::Map<Eigen::MatrixXd> (jac.data(),pos.size(),pos.size());
}
