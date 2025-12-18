// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisExpIm.h"
#include "basisSetDef.h"
#include "printOutput.h"
#include "useMatrix.h"
#include "basisSub.h"

BasisTrigon::BasisTrigon(const BasisSetDef &Def)
    :BasisIntegrable(Def.lowBound(),Def.upBound())
{
    if(lowBound()!=0. and abs(upBound()-2*math::pi)>1.e-12)
        ABORT(Str("trigonometric bases for now only on [0,2p], is: [","")+lowBound()+","+upBound()+"]");
    if(Def.funcs.find("expIm")!=0 and Def.funcs.find("cosSin")!=0 )DEVABORT("not a trigonometric definition: "+Def.funcs);
    _name=Def.funcs=="expIm" ? "expIm" : "CosSin"; //hmhm..rather inconsitent namings...

    std::vector<std::string> mstr;
    if(Def.funcs.find("[")!=std::string::npos)
        mstr=tools::splitString(tools::stringInBetween(Def.funcs,"[","]"),',');

    if(mstr.size()>1)
        for(std::string m: mstr)
            _mValues.push_back(tools::string_to_int(m));
    else {
        int mlow=0;
        if(mstr.size()==1)mlow=tools::string_to_int(mstr[0]);
        int curM=0;
        while(_mValues.size()<Def.order){
            if(_mValues.size()>0 or curM==mlow)_mValues.push_back(curM);
            // sorting will be 0,1,-1,2,-2,...
            curM=-curM;
            if(curM>=0)curM++;
        }
    }

    if(_mValues.size()<Def.order)
        ABORT(Sstr+"axis coefs="+Def.order+"> number of listed m's: "+Def.funcs);
    if(_mValues.size()>Def.order)
        PrintOutput::warning(Sstr+"axis coefs="+Def.order+"smaller than m's in list"+Def.funcs+" -- list will be truncated");

    _mValues.resize(Def.order);
}

BasisTrigon::BasisTrigon(std::string StrDefinition)
    :BasisIntegrable(0.,2.*math::pi)
{
    size_t colon=StrDefinition.find(":");
    _name=StrDefinition.substr(0,colon);
    //NOTE: for consistencys with old def, eventually switch to ExpIm etc.
    if(_name=="ExpIm")_name="expIm";
    else if(_name=="CosSin")_name="CosSin";
    else DEVABORT("string does not define ExpIm or CosSin, is: "+StrDefinition);
    std::vector<std::string> mStr=tools::splitString(StrDefinition.substr(colon+1),',');
    for(std::string m: mStr)_mValues.push_back(tools::string_to_int(m));
}

std::string BasisTrigon::str(int Level) const{
    if(size()<6)return Str(name(),"")+" {"+_mValues+"} ["+size()+"]";
    return Str(name(),"")+" {"+_mValues.front()+",...,"+_mValues.back()+"} ["+size()+"]";
}

bool BasisTrigon::operator==(const BasisAbstract& Other) const {
    const BasisTrigon* o=dynamic_cast<const BasisTrigon*>(&Other);
    if(o==0)return false;
    if(name()!=o->name())return false;
    return _mValues==o->_mValues;
}

std::string BasisTrigon::strDefinition() const {
    std::string s= dynamic_cast<const BasisExpIm*>(this) ? "ExpIm: " : "CosSin: ";
    for(int m: _mValues)s+=tools::str(m)+",";
    return s.substr(0,s.size()-1);
}

unsigned int BasisTrigon::order() const {
    unsigned int o=0;
    for(int k: _mValues)o=std::max(o,(unsigned int)std::abs(k));
    return 2*o+1;
}

void BasisTrigon::quadRule(int N, std::vector<double> &QuadX, std::vector<double> &QuadW) const{
    double dx=(upBound()-lowBound())/double(N);
    QuadX.clear();
    for (int k=0;k<N;k++)QuadX.push_back(lowBound()+dx*k);
    QuadW.assign(N,dx);
}

void BasisExpIm::valDer(const std::vector<std::complex<double> > &X,
                        std::vector<std::complex<double> > &Val,
                        std::vector<std::complex<double> > &Der, bool ZeroOutside) const{
    Val.clear();
    Der.clear();
    std::vector<std::complex<double> >expIx;
    for(std::complex<double> x: X)
        if(ZeroOutside and (x.real()<lowBound() or x.real()>upBound()))
            expIx.push_back(0.);
        else
            expIx.push_back(exp(std::complex<double>(0.,1.)*x));

    for(int m: _mValues)
        for(std::complex<double> y: expIx)
        {
            Val.push_back(std::pow(y,m)/sqrt(2.*math::pi));
            Der.push_back(std::complex<double>(0.,double(m))*Val.back());
        }
}
void BasisCosSin::valDer(const std::vector<std::complex<double> > &X,
                         std::vector<std::complex<double> > &Val,
                         std::vector<std::complex<double> > &Der, bool ZeroOutside) const{
    Val.clear();
    Der.clear();
    for(int m: _mValues){
        double qnrm= m==0 ? 1./sqrt(2.*math::pi) : 1./sqrt(math::pi) ;
        for(std::complex<double> x: X)
        {
            // not fast but clear
            if(m>0){
                Val.push_back(  sin(m*x.real())*qnrm);
                Der.push_back(m*cos(m*x.real())*qnrm);
            }
            else {
                Val.push_back(   cos(m*x.real())*qnrm);
                Der.push_back(-m*sin(m*x.real())*qnrm);
            }
        }
    }
}
