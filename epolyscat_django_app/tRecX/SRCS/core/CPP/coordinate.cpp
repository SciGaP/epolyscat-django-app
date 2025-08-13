// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coordinate.h"
#include "tools.h"
#include "constants.h"
#include "readInput.h"
#include "str.h"

#include <set>
using namespace physics;
using namespace math;
using namespace std;

map<string,Coordinate*> Coordinate::list;
void Coordinate::setUp(){
    // all available coordiante strings
    list["dum"]=  new Coordinate("dum" ,dum,false,false,J_one,    0.,DBL_MAX,"vector",false);
    list["BLANK"]=new Coordinate("NONE",dum,false,false,J_one,    0.,DBL_MAX,"vector",false);
    
    list["PEta" ]=new Coordinate("PEta",PEta, false, true, J_one, 0.,DBL_MAX,"sqrt[m%2{Phi}]*polynomial",false);
    list["PXi" ] =new Coordinate("PXi" ,PXi , false, true, J_one, 0.,DBL_MAX,"sqrt[m%2{Phi}]*polynomial",false);

    list["X"  ]=new Coordinate("X"   ,X,  true, true ,J_one,-DBL_MAX,DBL_MAX,"polynomial",false);
    list["Y"  ]=new Coordinate("Y"   ,Y,  true, true ,J_one,-DBL_MAX,DBL_MAX,"polynomial",false);
    list["Z"  ]=new Coordinate("Z"   ,Z,  true, true ,J_one,-DBL_MAX,DBL_MAX,"polynomial",false);
    list["Xh" ]=new Coordinate("Xh" , Xh, true, true ,J_one,       0,DBL_MAX,"polynomial",false);
    list["Phi"]=new Coordinate("Phi" ,Phi,false,false,J_one,      0.,  2.*pi,"cosSin",  true);
    list["The"]=new Coordinate("The" ,Th, false,false,J_one,   -pi/2,   pi/2,"cosSin",  true);
    list["CTh"]=new Coordinate("CTh" ,CTh,false,false,J_one,     -1,       1,"polynomial",false);
    list["Xi" ]=new Coordinate("Xi"  ,Xi, false,false,J_one,      1.,DBL_MAX,"polynomial",false);
    list["Eta"]=new Coordinate("Eta" ,Eta,false,false,J_one,     -1.,     1.,"assocLegendre{Phi}",true);
    list["L"  ]=new Coordinate("L",   L,  false,false,J_one,     -1.,     1.,"polynomial",false);
    list["M"  ]=new Coordinate("M",   M,  false,false,J_one,      0.,DBL_MAX,"vector",false);
    list["Rho"]=new Coordinate("Rho" ,Rho,false,true ,J_val,      0.,DBL_MAX,"polynomial",false);
    list["Rhn"]=new Coordinate("Rhn" ,Rhn,true, true ,J_one,      0.,DBL_MAX,"sqrt*polynomial",false);
    list["R"  ]=new Coordinate("R"   ,R,  false,true ,J_squ,      0.,DBL_MAX,"polynomial",false);
    list["Rn" ]=new Coordinate("Rn"  ,Rn, true, true, J_one,      0.,DBL_MAX,"polynomial",false);
    list["K"  ]=new Coordinate("K"   ,K,  false,false,J_one,     -1.,     1.,"legendre",false);
    list["Idx"]=new Coordinate("Idx" ,Idx,false,false,J_one,      0.,DBL_MAX,"vector",false);
    list["Vec"]=new Coordinate("Vec" ,Vec,false,false,J_one,      0.,DBL_MAX,"vector",false);
    list["kGrid"]=new Coordinate("kGrid",Vec,false,false,J_one,-DBL_MAX,DBL_MAX,"vector",false);

    list["Ion"]=new Coordinate("Ion" ,Ion,false,false,J_one,      0.,DBL_MAX,"CIion",false);
    list["Neutral"]=new Coordinate("Neutral" ,Neut,false,false,J_one,0.,DBL_MAX,"CIneut",false);
    list["Hybrid"]=new Coordinate("Hybrid",Vec,false,false,J_one,      0.,DBL_MAX,"vector",false);;
    list["Orbital"]=new Coordinate("Orbital" ,Orbital,false,false,J_one,0.,DBL_MAX,"NO_VALUE",false);
    list["Ndim"]  =new Coordinate("Ndim" ,Vec,false,false,J_one,0.,DBL_MAX,"NO_VALUE",false);

    // various aliases
    list["specRn"]=list["dum"];
    list["specX"]=list["dum"];
    list["specY"]=list["dum"];
    list["specZ"]=list["dum"];
    list["kR"]=list["dum"];
    list["kX"]=list["dum"];
    list["kY"]=list["dum"];
    list["kZ"]=list["dum"];
    list["v/d"]=list["dum"];
    list["nSurface"]=list["dum"];
    list["kRn"]=list["kGrid"];
    list["Channel"]=list["Vec"];
}
void Coordinate::cleanUp(){
    vector<string> alias({"v/d","nSurface","Channel"});
    for (auto p: list){
        if(find(alias.begin(),alias.end(),p.first)==alias.end()){
            if(p.first.find("k")!=0 and p.first.find("spec")!=0){
                delete p.second;
                p.second=0;
            }
        }
    }
}

Coordinate Coordinate::fromString(std::string String0){
    return Coordinate(String0);
}

string Coordinate::kind(string Name){return Coordinate(Name).cString;}

Coordinate::Coordinate(std::string String){
    // strip trailing numbers from string
    std::string string=String.substr(0,String.find("periodic"));
    string=string.substr(0,String.find_first_of("0123456789"));
    if(list.size()==0)Coordinate::setUp();
    if(list.count(string)==0)ABORT("coordinate not defined: "+string);
    operator=(*list[string]);
    if(list[string]->cString=="undefined")ABORT("coordinate not defined: "+string);
    if(String.find("periodic")!=std::string::npos) {
        if(string!="X" and string!="Y" and string!="Z")
            ABORT("only X,Y,Z can be defined periodic, not "+String);
        _periodic=true;
    }
}

std::string Coordinate::name() const { return cString; }
std::string Coordinate::defaultFunction() const {return function;}

int Coordinate::automaticOrder(std::string String, int NCoef){
        if(String=="dum")return 0;
        string string0=String.substr(0,String.find_first_of("0123456789"));
        vector<string> orderAll={"Channel","Vec","Orbital","Ion","Neutral","Hybrid","CTh","Phi","Eta","Ndim"};
        if(std::find(orderAll.begin(),orderAll.end(),string0)!=orderAll.end())return NCoef;
        return -1;
}

static std::vector<std::string> discreteBases({"vector","CIon","CIneut","NO_VALUE"});
bool Coordinate::isDiscrete(std::string Name){
    if(list.count(Name)==0)return false; // numbered axes are continuous as a rule ??

    return std::find(discreteBases.begin(),discreteBases.end(),list[Name]->function)!=discreteBases.end();
}
