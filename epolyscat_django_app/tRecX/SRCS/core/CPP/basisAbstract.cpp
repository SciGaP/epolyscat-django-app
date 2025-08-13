// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisAbstract.h"

#include "readInput.h"
#include "printOutput.h"

#include "basisSetDef.h"
#include "basisAbstract.h"
#include "basisNdim.h"
#include "basisDvr.h"
#include "basisDvrSqrt.h"
#include "basisSqrtDVR.h"
#include "basisGridQuad.h"
#include "basisAbstract.h"
#include "basisBesselCoulomb.h"
#include "basisEigen.h"
#include "basisExpIm.h"
#include "basisIntegrableSqrt.h"
#include "basisIntegrableRpow.h"
#include "basisVector.h"
#include "basisOrbitalNumerical.h"
#include "basisChannel.h"
#include "basisSub.h"
#include "threads.h"
#include "basisAssocLeg.h"
#ifdef _USE_HACC_
#include "basisfunctioncineutral.h"
#include "basisfunctionciion.h"
#include "basisCI.h"
#endif

#include "basisIntegrable.h"
#include "basisSub.h"
#include "basisGrid.h"
#include "basisGridQuad.h"
#include "basisNdim.h"
#include "basisHybrid.h"
#include "basisMO.h"

static std::map<std::string,const BasisAbstract*> _listOfBases;

const BasisNdim* BasisAbstract::ndim() const {return dynamic_cast<const BasisNdim*>(this);}
const BasisGridQuad* BasisAbstract::gridQuad() const{return dynamic_cast<const BasisGridQuad*>(this);}
const BasisSub* BasisAbstract::sub() const{return dynamic_cast<const BasisSub*>(this);}
const BasisOrbital* BasisAbstract::orbital() const{return dynamic_cast<const BasisOrbital*>(this);}
const BasisHybrid* BasisAbstract::hybrid() const{return dynamic_cast<const BasisHybrid*>(this);}
const BasisIntegrable* BasisAbstract::integrable() const {
  const BasisIntegrable* b=dynamic_cast<const BasisIntegrable*>(this);
    if(b==0)b=dynamic_cast<const BasisDVR*>(this);
    if(b==0)return 0;
    if(b->isDVR())return b;
    if(b->isGrid())return 0;
    if(b->isIndex())return 0;
    return b;
}
const BasisGrid* BasisAbstract::grid() const {
    const BasisGrid* b=dynamic_cast<const BasisGrid*>(this);
    if(b)return b;
    if(not BasisSub::superBas(this)->isGrid())return 0;
    return BasisGrid::factory(this);
}

#ifdef _USE_HACC_
const BasisCI* BasisAbstract::ci() const{return dynamic_cast<const BasisCI*>(this);}
#endif

const BasisAbstract* BasisAbstract::factory(const BasisSetDef & Def){

    std::string def;
    const BasisAbstract* bas(0);
    if(0!=(bas=BasisNdim::factory(Def.funcs)));
    else if(Def.funcs=="polynomial")bas=new BasisDVR(Def);
    else if(Def.funcs.find("DVR: jacobi[")==0)bas=new BasisDVR(Def);
    else if(Def.funcs.find("polExp")==0)bas=new BasisDVR(Def);
    else if(Def.funcs.find("besselCoulomb")==0)bas=new BasisBesselCoulomb(Def);
    else if(Def.funcs.find("Eigenbasis")==0)bas=new BasisEigen(Def);
    else if(Def.funcs.find("expIm")==0)bas=new BasisExpIm(Def);
    else if(Def.funcs.find("cosSin")==0)bas=new BasisCosSin(Def);
    else if(Def.funcs.find("vector")==0)bas=new BasisVector(Def.order);
    else if(Def.funcs.find("Orbital")==0)bas=new BasisOrbitalNumerical(Def);
    else if(Def.funcs.find("Channel")==0)bas=new BasisChannel(Def.funcs,Def.order,int(Def.lowBound()));
    else if(Def.name().find("assocLegendre")==0)bas=new BasisAssocLeg(BasisAssocLeg::strDefinition(Def));
#ifdef _USE_HACC_
    else if(Def.funcs.find("CIneut")==0)bas=new BasisFunctionCINeutral(BasisFunctionCINeutral::strDefinition(Def));
    else if(Def.funcs.find("CIion")==0)bas=new BasisFunctionCIion(BasisFunctionCIion::strDefinition(Def));
    else if(Def.funcs.find("CI[")==0)bas=new BasisCI(BasisCI::strDefinition(Def));
#endif

    else if(""!=(def=BasisIntegrableSqrt::inputToStrDefinition(Def)))return factory(def);
    else if(""!=(def=BasisIntegrableRpow::inputToStrDefinition(Def)))return factory(def);

    if(not bas)bas=BasisEigen::factory(Def.funcs+":"+tools::str(int(Def.lowBound()))+":"+tools::str(int(Def.lowBound()+Def.order)));

    // use old style BasisSet
    if(bas==nullptr)DEVABORT("cannot interprete basis intput: "+Def.str());
    def=bas->strDefinition();

    // PolarOff can only be read from input
    if(def.find("PolarOff")!=std::string::npos)return bas;
    if(Def.funcs.find("Channel")!=std::string::npos)return bas;

    delete bas;
    return factory(def);
}

std::vector<const BasisAbstract*> BasisAbstract::select(std::function<bool (const BasisAbstract*)> Criterion){
    std::vector<const BasisAbstract*> res;
    for(auto p: _listOfBases)
        if(Criterion(p.second))
            res.push_back(p.second);
    return res;
}

const BasisAbstract* BasisAbstract::factory(const std::string & Def){
    //NOTE: the specialized factories should all be eliminated
    std::string def=tools::cropString(Def);
    if(_listOfBases.count(def))return _listOfBases[def];
    const BasisAbstract* bas(0);
    if(0==def.find("Thread:"))      bas=new BasisThread(def);
    if(0==def.find("BasisSub:"))    bas=new BasisSub(def);
    if(0==def.find("ExpIm:"))       bas=new BasisExpIm(def);
    if(0==def.find("CosSin:"))      bas=new BasisCosSin(def);
    if(0==def.find("AssociatedLegendre:"))bas=new BasisAssocLeg(def);
    if(0==def.find("Grid:"))        bas=BasisGrid::factory(def);
    if(0==def.find("GridQuad:"))    bas=BasisGridQuad::factory(def);
    if(0==def.find("DVR:"))         bas=new BasisDVR(def);
    if(0==def.find("NONE:"))        bas=BasisVector::factory(def);
    if(0==def.find("Coef:"))        bas=BasisVector::factory("Vector:"+def.substr(5));
    if(0==def.find("Hybrid:"))      bas=new BasisHybrid(tools::string_to_int(def.substr(7)));
    if(0==def.find("Vector:"))      bas=BasisVector::factory(def);
    if(0==def.find("Orbital"))      bas=new BasisOrbitalNumerical(def);
    if(0==def.find("Rpow"))         bas=new BasisIntegrableRpow(def);
#ifdef _USE_HACC_
    if(0==def.find("CIneutral:"))   bas=new BasisFunctionCINeutral(def);
    if(0==def.find("CIion:"))       bas=new BasisFunctionCIion(def);
    if(0==def.find("CI["))          bas=new BasisCI(def);
#endif
    if(0==def.find("MO[")){
        bas=dynamic_cast<const BasisMO*>(VectorValuedFunction::get(tools::stringInBetween(def,"[","]",true)));
        if(not bas)DEVABORT("failed to retrieve "+def);
    }

    if(0==def.find("sqrt")){
        if(def.find("DVR:")!=std::string::npos)bas=new BasisSqrtDVR(def);
        else                                   bas=new BasisIntegrableSqrt(def);
    }

    if(0==bas)bas=BasisEigen::factory(def);

    if(0==bas)DEVABORT("construction from StrDefinition not implemented: "+def.substr(0,80)+"...");

    if(dynamic_cast<const BasisIntegrable*>(bas) and
            ReadInput::main.flag("DEBUGplotBasis","plot all basis functions into run output directory"))
        dynamic_cast<const BasisIntegrable*>(bas)->plot(ReadInput::main.output());

    return _listOfBases[def]=bas;
}
bool BasisAbstract::isAbsorptive() const {return integrable()?integrable()->isAbsorptive():false;}

const BasisAbstract * BasisAbstract::remove(const std::vector<int> &RemoveK) const{
    std::vector<int> subset;
    for(size_t k=0;k<size();k++)
        if(std::find(RemoveK.begin(),RemoveK.end(),k)==RemoveK.end())
            subset.push_back(k);
    return BasisAbstract::factory(BasisSub::strDefinition(this,subset));
}

bool BasisAbstract::operator==(const BasisAbstract &Other) const{
    if(this==&Other)return true;
    if(size()!=Other.size())return false;
    if(isGrid())return BasisGrid::factory(this)->operator==(Other);
    if(dynamic_cast<const BasisTrigon*>(this))return dynamic_cast<const BasisTrigon*>(this)->operator==(Other);
    if(dynamic_cast<const BasisOrbitalNumerical*>(this)!=0)return dynamic_cast<const BasisOrbitalNumerical*>(this)->operator==(Other);
    if(_name=="NONE" and Other._name=="NONE")return true;
    if((dynamic_cast<const BasisThread*>(this)==0)!=(dynamic_cast<const BasisThread*>(&Other)==0))return false;
    if(dynamic_cast<const BasisSub*>(this))return dynamic_cast<const BasisSub*>(this)->operator==(Other);
    if(dynamic_cast<const BasisVector*>(this))return dynamic_cast<const BasisVector*>(this)->operator==(Other);

    // fall back to comparing strDefinition (may be slow, implement direct comparison recommended)
    if(BasisAbstract::strDefinition()!=strDefinition() and Other.BasisAbstract::strDefinition()!=Other.strDefinition())
        return strDefinition()==Other.strDefinition();
    DEVABORT(Sstr+"comparison not implemented:"+str()+"??"+Other.str());
}

std::string BasisAbstract::str(int Level) const{
    if(dynamic_cast<const BasisNdim*>(this)!=0)return dynamic_cast<const BasisNdim*>(this)->str(Level);
    if(dynamic_cast<const BasisDVR*>(this)!=0)return dynamic_cast<const BasisDVR*>(this)->str(Level);
    return Str(_name," ")+size();
}
std::string BasisAbstract::strDefinition() const{
    if(name()=="NONE")return "NONE:";
    return Str(name(),"")+","+size();
}
