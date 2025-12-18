// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "parameters.h"
#include "readInput.h"
#include "algebra.h"
#include "constants.h"

using namespace std;

deque<Parameters> Parameters::table;
deque<Updatable*> Parameters::updatables;

double Parameters::lastUpdateTime=-DBL_MAX;
Parameters Parameters::noParameter;
bool Parameters::parsForcedTo1=false;

Updatable::~Updatable() {}

std::complex<double> * Parameters::pointer(std::string Name){
    if(Name=="")return &(tableParameter(Name).plusValue);
    if(Name[0]=='-'){return &(tableParameter(Name.substr(1)).minusValue);};
    if(Name[0]=='+'){return &(tableParameter(Name.substr(1)).plusValue );};
    return &(tableParameter(Name).plusValue);
}
void Parameters::defaults(){
    for(int k=0;k<table.size();k++)
        if(table[k].name=="")return;
    table.push_back(Parameters("",1.,0,false)); // empty string is added by force
    add("i",std::complex<double>(0.,1),0);
}

/// parFunction pointer
void Parameters::add(string Name, const complex<double> CurrentValue, parFunction Update){
    if(Name=="")ABORT("must not add empty name - is set by default");
    if(table.size()==0)defaults();
    if(Name[0]=='+' or Name[0]=='-')ABORT("do not include sign in parameter name "+Name);
    bool exists=false;
    // search in standard parameter list
    for(unsigned int n=0;n<table.size();n++){
        if(table[n].name==Name){
            exists=true;
            if(table[n].updateFunction!=Update) ABORT("conflicting redefinition of parameter function: "+Name);
            if(table[n].updateFunctionArg!=0) ABORT("conflicting redefinition of parameter function: "+Name);
            if(table[n].algebra!=0) ABORT("conflicting redefinition of parameter function: "+Name);
            if((Update==0) and (table[n].plusValue!=CurrentValue))ABORT("conflicting redefinition of parameter value: "+Name);
        }
    }
    complex<double>val=CurrentValue;
    if(Update!=0)val=1.;
    if(not exists){
        table.push_back(Parameters(Name,val,Update,false));
    }
}
/// parFunction pointer
void Parameters::add(const Algebra* Alg){
    std::string Name=Alg->definition();
    if(Name=="")ABORT("must not add empty name - is set by default");
    if(table.size()==0)defaults();
    if(Name[0]=='+' or Name[0]=='-')ABORT("do not include sign in parameter name "+Name);
    bool exists=false;
    // search in standard parameter list
    for(unsigned int n=0;n<table.size();n++){
        if(table[n].name==Name){
            exists=true;
            if(not table[n].algebra or
                    table[n].algebra->definition()!=Alg->definition())
                ABORT("conflicting redefinition of parameter function: "+Name);
        }
    }
    if(not exists){
        complex<double>val=Alg->val(0.);
        table.push_back(Parameters());
        table.back().name=Alg->definition();
        table.back().algebra=Alg;
        table.back().plusValue=val;
        table.back().minusValue=-val;
    }
}

void Parameters::setSpecial(){

    if(table.size()==0)defaults();

    // polar angles and X,Y,Z
    if(Parameters::isDefined("radius") and Parameters::isDefined("theta"))
    {
        Parameters::addResetable("Z[radius,theta]");
        if(Parameters::isDefined("phi")){
        }
    }

    // lambda and intensity and Omega,A,F
    if(Parameters::isDefined("lambda"))Parameters::addResetable("Omega[lambda]");
    if(Parameters::isDefined("intensity"))Parameters::addResetable("F[intensity]");
    if(Parameters::isDefined("intensity") and Parameters::isDefined("lambda"))Parameters::addResetable("A[lambda,intensity]");
}

void Parameters::updateSpecial(){
    // polar angles and X,Y,Z
    double the=DBL_MAX,rad=DBL_MAX,phi=DBL_MAX;
    if(Parameters::isDefined("radius"))rad=Parameters::pointer("radius")->real();
    if(Parameters::isDefined("theta")) the=Parameters::pointer("theta")->real();
    if(Parameters::isDefined("phi")) phi=Parameters::pointer("phi")->real();
    if(rad!=DBL_MAX and the!=DBL_MAX){
        if(not Parameters::isDefined("Z[radius,theta]"))Parameters::addResetable("Z[radius,theta]");
        Parameters::reset("Z[radius,theta]",rad*cos(the));
        if(phi!=DBL_MAX){
            Parameters::reset("X[radius,phi,theta]",rad*sin(the)*cos(phi));
            Parameters::reset("Y[radius,phi,theta]",rad*sin(the)*sin(phi));
        }
    }

    // lambda and intensity and Omega,A,F
    double lam=-1.,inten=-1.;
    double omega;
    if(Parameters::isDefined("lambda"))lam=Parameters::pointer("lambda")->real();
    if(Parameters::isDefined("intensity"))inten=Parameters::pointer("intensity")->real();
    if(lam>0.){
        omega=2*math::pi/(physics::a_finestructure*lam);
        Parameters::reset("Omega[lambda]",omega);
        if(inten>=0)Parameters::reset("A[lambda,intensity]",sqrt(inten)/omega);
    }
    if(inten>=0)Parameters::reset("F[intensity]",sqrt(inten));
}

/// parFunction pointer
void Parameters::addResetable(string Name, const complex<double> CurrentValue){
    if(table.size()==0)defaults();
    if(Name[0]=='+' or Name[0]=='-')ABORT("do not include sign in parameter name "+Name);
    bool exists=false;
    // search in standard parameter list
    for(unsigned int n=0;n<table.size();n++){
        if(table[n].name==Name){
            exists=true;
            if(table[n].updateFunction!=0) ABORT("conflicting redefinition of parameter function: "+Name);
            if(table[n].updateFunctionArg!=0) ABORT("conflicting redefinition of parameter function: "+Name);
            if(table[n].algebra!=0) ABORT("conflicting redefinition of parameter function: "+Name);
            if(not table[n].resetable)ABORT("conflicting redefinition of parameter value: "+Name);
        }
    }
    complex<double>val=CurrentValue;
    if(not exists)table.push_back(Parameters(Name,val,0,true));
}

bool Parameters::isDefined(std::string Name){
    for(auto p: table)
        if(p.name==tools::cropString(Name))return true;
    return false;
}

bool Parameters::isFunction(std::string Name){
    if(Name[0]=='-' or Name[0]=='+')
        return tableParameter(Name.substr(1)).updateFunction!=0
                or tableParameter(Name.substr(1)).algebra!=0
                or tableParameter(Name.substr(1)).resetable;
    return tableParameter(Name).updateFunction!=0 or tableParameter(Name).algebra!=0
            or tableParameter(Name).updateFunctionArg!=0 or tableParameter(Name).resetable ;
}

void Parameters::update(){
    for(Parameters &p: table){
        if(p.resetable){
            Algebra upd(p.name);
            if(upd.isAlgebraOfUpdatables()){
                p.plusValue=upd.val(0.);
                p.minusValue=-p.plusValue;
            }
        }
    }
}

void Parameters::show(){
    std::cout<<"Available parameters:\nName = currentValue"<<std::endl;
    for (unsigned int k=0;k<table.size();k++){
        std::cout<<"'"<<table[k].name<<"' = "<<table[k].plusValue;
        if(table[k].updateFunction!=0 or table[k].updateFunctionArg!=0 or table[k].algebra!=0)std::cout<<" [Function(t)]";
        else if(Algebra(table[k].name).isAlgebraOfUpdatables())std::cout<<" (updatabel)";
        else if(table[k].resetable)std::cout<<" (resetable)";
        std::cout<<std::endl;
    }
}

void Parameters::reset(string Name, const complex<double> CurrentValue){
    if(Name[0]=='+' or Name[0]=='-')ABORT("do not include sign in parameter name "+Name);
    // search in standard parameter list
    for(unsigned int n=0;n<table.size();n++){
        if(table[n].name==Name){
            if(not table[n].resetable)ABORT("for resetable parameter, add "+Name+" by \"addResetable(...)\" (performance penalty!)");
            table[n].plusValue=CurrentValue;
            table[n].minusValue=CurrentValue;
            update(); // parameters that may be dependent on this one (updatable algebras)
            return;
        }
    }
    show();
    ABORT("parameter not in table: "+Name+", call Parameters::add(...) first");
}

/// FunctionOneArg pointer
void Parameters::add(string Name){
    if(table.size()==0)defaults();
    if(Name[0]=='+' or Name[0]=='-')ABORT("do not include sign in parameter name "+Name);
    // search in standard parameter list
    for(unsigned int n=0;n<table.size();n++){
        if(table[n].name==Name)return;
        table.push_back(Parameters(Name));
    }
}

Parameters & Parameters::tableParameter(const string Name,bool Abort)
{
    if(table.size()==0){
        defaults();
        setSpecial();
    }

    for(unsigned int n=0;n<Parameters::table.size();n++)
        if(Parameters::table[n].name==tools::cropString(Name))return Parameters::table[n];
    if(FunctionOneArg::get(Name,false)!=0){
        // try adding from general functions
        Parameters::table.push_back(Parameters(Name));
    } else if (Name.find_first_not_of("+-.0123456789e")>=Name.length() and tools::string_to_double(Name)!=0.){
        Parameters::add(Name,tools::string_to_double(Name));
    } else if (Algebra(Name).isAlgebraOfUpdatables()){
        Parameters::addResetable(Name,Algebra(Name).val(0.));
    } else if (Algebra(Name).isAlgebraOfConsts()){
        Parameters::add(Name,Algebra(Name).val(0.));
    } else if (Algebra(Name).isAlgebra()){
        Algebra * a=new Algebra(Name);
        Parameters::add(a);
    }
    Algebra::failures=""; // failure here is admitted
    if(table.back().name==Name)return table.back();

    if(Abort){
        show();
        ABORT("no operator parameter defined for string: \""+Name
              +"\"\nadd by Parameters::add(string,value,functionPointer) "+
              +"\n    or Parameters::add(string) "
              +"\n    or Parameters::setSpecial() for dependent parameters"
              +"\n    or Algebra::specialConstants"
              +"\n    or correctly formed algebraic expression of the above"
              );
    } else {
        return noParameter;
    }
}

void Parameters::read(string Cat, string Name, string Default, ReadInput &in){
    double val;
    in.read(Cat,Name,val,Default,"add value as "+Cat+":"+Name+" to parameters");
    add(Cat+":"+Name,val,0);
}

void Parameters::add(Updatable* ptUpdatable)
{
    for (deque<Updatable*>::iterator it(updatables.begin()); it!=updatables.end(); ++it) {
        if (*it==ptUpdatable) {return;}
    }
    updatables.push_back(ptUpdatable);
}

TIMER(updatePar,)
TIMER(updatable,)
void Parameters::update(const double time, bool Special) {

    if(parsForcedTo1)DEVABORT("Parameters::updateToOne() was not followed by restoreToTime()");

    lastUpdateTime=time;
    for(unsigned int n=0;n<table.size();n++){
        if     (table[n].algebra)          table[n].plusValue=table[n].algebra->val(time).real();
        else if(table[n].updateFunctionArg)table[n].plusValue=table[n].updateFunctionArg->val(time);
        else if(table[n].updateFunction)   table[n].plusValue=(table[n].updateFunction)(time);
        table[n].minusValue=-table[n].plusValue;
    }
    for (deque<Updatable*>::iterator it(updatables.begin()); it!=updatables.end(); ++it) {
        (*it)->update(time);
    }
    if(Special)updateSpecial();
}


void Parameters::updateToOne() {
    parsForcedTo1=true;
    for(unsigned int n=0;n<table.size();n++)
        if(table[n].algebra!=0 or table[n].updateFunction!=0 or table[n].updateFunctionArg!=0)
        {
            table[n].plusValue =1.;
            table[n].minusValue=1.;
        }
}




