// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "units.h"
#include "constants.h"
#include "tools.h"
#include "printOutput.h"

using namespace std;
using namespace math;
using namespace physics;

string Units::sep="|";
map<string,double>Units::uniTab;
map<string,string>Units::systems,Units::aka;

void Units::setDefault(string Def){
    if(not tools::hasKey(systems,string("SI")))standardUnits();
    if(not tools::hasKey(systems,Def))
        ABORT("unit system not found: "+Def+"\n"+tools::listMapKeys(systems,", "));
    aka["DEFAULT_SYSTEM"]=Def;
}

void Units::standardUnits(){
    if(tools::hasKey(systems,string("standardUnits")))return;
    systems["standardUnits"]="set";
    //    if(not tools::hasKey(systems,string("DEFAULT_SYSTEM")))aka["DEFAULT_SYSTEM"]="SI";

    // systems
    addUnitSystem("SI",1.,1.,1.,1.,4*pi*1.e-7,"m|length,s|time,kg|mass,C|charge,A|current,V/m|eField");
    addUnitSystem("atomic units"+sep+"au",bohr_radius,electron_mass,bohr_radius/(speed_of_light*a_finestructure),
                  proton_charge,4.*pi*pow(a_finestructure,2),"length|Bohr");
    addUnitSystem("relativistic"+sep+"rel",h_bar/(1.e6*proton_charge),1.e6*proton_charge,
                  h_bar/(1.e6*proton_charge*speed_of_light),sqrt(4*pi*a_finestructure),0.5/pi);
    addUnitSystem("cm-g-sec(ESU)"+sep+"cgs",1e-2,1e-3,1.,1e-1/speed_of_light,4*pi/pow(speed_of_light*1.e4,2));

    // individual units
    addUnit("Ry"+sep+"energy",0.5*uniTab["au"+sep+"energy"]);
    addUnit("eV"+sep+"energy",proton_charge);
    addUnit("W/cm2"+sep+"intensity",1.e4);
    addUnit("nm"+sep+"length",1.e-9);
}

void Units::addUnitSystem(string System, double Length, double Mass, double Time, double Charge, double Mu0, string Names){

    string name=tools::cropString(System);
    size_t col=name.find(sep);
    if(col==string::npos){
        col=name.length();
        name=name+sep+name;
    }
    string longName=name.substr(0,col-1),shortName=name.substr(col+1);
    if(longName=="" or shortName=="")ABORT("specify non-empty long and short unit names, format 'long:short', input is "+System);

    systems[shortName]=longName;
    aka[shortName]=shortName;

    addUnit(shortName+sep+"length",Length);
    addUnit(shortName+sep+"mass",  Mass);
    addUnit(shortName+sep+"time",  Time);
    addUnit(shortName+sep+"charge",Charge);
    addUnit(shortName+sep+"mu0",Mu0);

    // derived units
    addUnit(shortName+sep+"velocity",uniTab[shortName+sep+"length"]/uniTab[shortName+sep+"time"]);
    addUnit(shortName+sep+"energy",pow(uniTab[shortName+sep+"velocity"],2)*uniTab[shortName+sep+"mass"]);
    addUnit(shortName+sep+"frequency",1./uniTab[shortName+sep+"time"]);
    addUnit(shortName+sep+"ep0",pow(uniTab[shortName+sep+"velocity"],2)/(pow(speed_of_light,2)*uniTab[shortName+sep+"mu0"]));
    addUnit(shortName+sep+"action",uniTab[shortName+sep+"energy"]*uniTab[shortName+sep+"time"]);

    if(shortName=="SI"){
        addUnit(shortName+sep+"eField",1.);
        addUnit(shortName+sep+"intensity",1.);
    } else {
        addUnit(shortName+sep+"eField",uniTab[shortName+sep+"charge"]/pow(uniTab[shortName+sep+"length"],2)
                /(4*pi*uniTab["SI"+sep+"ep0"]));
        addUnit(shortName+sep+"intensity",pow(uniTab[shortName+sep+"eField"],2)/2*speed_of_light*uniTab["SI"+sep+"ep0"]);
    }

    // add extra names to aka-table
    string names=Names;
    tools::cropString(names);
    while (names!=""){
        size_t colon=names.find(sep),end=min(names.find(","),names.length()-1);
        if(colon==0 or colon==string::npos or end==colon)ABORT("corrupted list of names: "+Names);
        aka[names.substr(0,colon)]=shortName+sep+names.substr(colon+1,end-colon-1);
        names=names.substr(end+1);
    }
}

void Units::addUnit(string NameDimension, double Value){
    standardUnits();

    if(NameDimension.find(sep)==string::npos)
        ABORT("need format unitsName"+string(sep)+"dimensionKind, is:"+NameDimension);
    if(tools::hasKey(uniTab,NameDimension))
        if(uniTab[NameDimension]!=Value)
            ABORT("conflicting definition of unit: present "+tools::str(Value,14)+" vs. previous "+tools::str(uniTab[NameDimension],14));
    string name=NameDimension.substr(0,NameDimension.find(sep));
    aka[NameDimension]=NameDimension; // add full name

    // if first part is not system name, add first part as synonym
    if(not tools::hasKey(systems,name)){
        aka[NameDimension.substr(0,NameDimension.find(sep))]=NameDimension;
        //        cout<<"non-standard AKA: "+NameDimension.substr(0,NameDimension.find(sep))+"///"+aka[NameDimension.substr(0,NameDimension.find(sep))]<<endl;
    }
    uniTab[NameDimension]=Value;
}

double Units::convert(double Val, string InUnits, string OutUnits){
    return convert(vector<double>(1,Val),InUnits,OutUnits)[0];
}

bool Units::isDefined(string Name){
    return aka.count(Name)==1;
}
vector<complex<double> > Units::convert(vector<complex<double> > Val,  string InUnits, string OutUnits){
    vector<complex<double> > d;
    for(complex<double> c: Val)d.push_back(complex<double>(
                                               convert(c.real(),InUnits,OutUnits),
                                               convert(c.imag(),InUnits,OutUnits))
                                           );
    return d;
}
vector<double> Units::convert(vector<double> Val,  string InUnits, string OutUnits){
    // no need to convert
    InUnits=tools::cropString(InUnits);
    OutUnits=tools::cropString(OutUnits);
    if(InUnits==OutUnits)return Val;

    standardUnits();

    // find units in aka list
    string inp=aka[InUnits],out=aka[OutUnits];

    string inpDim=inp.substr(min(inp.find(sep),inp.length()));
    string outDim=out.substr(min(out.find(sep),out.length()));

    // empty strings may indicate unit system
    if(inpDim==""){
        inpDim=outDim;
        inp+=inpDim;
    }
    if(outDim==""){
        outDim=inpDim;
        out+=outDim;
    }
    if(inpDim=="")ABORT("cannot infer dimension from input or output: "+InUnits+" -> "+OutUnits);
    string errMess;
    if(not tools::hasKey(uniTab,inp))errMess="\nERROR: Input units not defined: "+InUnits;
    if(not tools::hasKey(uniTab,out))errMess="\nERROR: Output units not defined: "+OutUnits;
    if(errMess!=""){
        if(errMess.find("DEFAULT_SYSTEM"))errMess+="\n       set a default unit system by Units::setDefault(system-name)";
        errMess+="\n       select from list above";
        ABORT("available units:\nsystem-name|dimension\n---------------\n"+tools::listMapKeys(uniTab,"\n")+errMess);
    }

    if(uniTab.size()==0)standardUnits();
    if(inpDim!=outDim)ABORT("cannot convert "+InUnits+" -> "+OutUnits+" (aka: "+inp+" -> "+out+")");

    vector<double>conVal;
    for(unsigned int k=0;k<Val.size();k++)
        conVal.push_back(Val[k]*uniTab[inp]/uniTab[out]);
    return conVal;
}

void Units::test(){
    cout<<tools::listMap(uniTab);
    setDefault("SI");
    check("au"+string(sep)+"eField",5.14220631797628e11,11);
    check("au"+string(sep)+"intensity",3.509e20,3);
    check("au"+string(sep)+"energy",27.211*convert(1,"eV"+string(sep)+"energy"),4);
    check("Ry",convert(1.,"au"+string(sep)+"energy")*0.5,14);
    check("m",1,14);
    check("cgs"+string(sep)+"eField",speed_of_light*1.e-4,14);
    check("cgs"+string(sep)+"charge",1.e-1/speed_of_light,14);

    setDefault("au");
    check("Ry",0.5,14);

    setDefault("SI");
    check("cgs"+string(sep)+"charge",1.e-1/speed_of_light,14);
    PrintOutput::message("Units::test passed");
}
void Units::check(const string Unit, double Val, int Digits){
    double conVal=convert(1.,Unit);
    if(abs(conVal/Val-1)>pow(0.1,Digits))
        ABORT(Unit+" incorrect by > "+tools::str(Digits)+" Digits, (value,ratio)=("+tools::str(conVal)+", "+tools::str(conVal/Val,10)+")");
    cout<<"OK units: "+Unit<<endl;
}




