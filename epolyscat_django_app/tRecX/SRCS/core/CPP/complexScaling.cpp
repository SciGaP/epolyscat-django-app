// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "complexScaling.h"
#include "tools.h"
#include "basisMat.h"
#include "readInput.h"

using namespace std;

vector<string> ComplexScaling::names;

/// @cond DEV

class pmlEta:public basisMatFunc{
    friend class ComplexScaling;
    friend class pmlRexpSig;
    friend class pmlRexpRoot;
    ComplexScaling comSca;
    pmlEta(const ComplexScaling & ComSca):basisMatFunc("JpmlEta"),comSca(ComSca){}
    double r0lower() const {return comSca._r0lower;}
    double r0upper() const {return comSca._r0upper;}
public:
    complex<double> operator()(complex<double> z) const {
        if(z.real()<comSca._r0lower or z.real()>comSca._r0upper)return comSca.eta;
        else return 0;
    }
    bool operator==(const basisMatFunc & Other) const {
        if(name!=Other.name)return false;
        /// note: dynamic_cast may require compiler option on some compilers
        const pmlEta * o=dynamic_cast<const pmlEta*>(&Other);
        if(o==0)ABORT("two different basisMatFunc instances with same name: "+o->name);
        return o->comSca._r0lower==comSca._r0lower and o->comSca._r0upper==comSca._r0upper and o->comSca.eta==comSca.eta;}
};
class pmlEtaRoot:public basisMatFunc{
    friend class ComplexScaling;
    pmlEta eta;
    pmlEtaRoot(const pmlEta & Eta):basisMatFunc("JpmlEtaRoot"),eta(Eta){}
public:
    complex<double> operator()(complex<double> z) const {return 0.25*eta(z)/sqrt(z);}
    bool operator==(const basisMatFunc & Other) const {
        if(name!=Other.name)return false;
        const pmlEtaRoot * o=dynamic_cast<const pmlEtaRoot*>(&Other);
        if(o==0)ABORT("two different basisMatFunc instances with same name: "+o->name);
        return o->eta==eta;}
};
class pmlRexpSig:public basisMatFunc{
    friend class ComplexScaling;
    pmlEta eta;
    bool plus;
    pmlRexpSig(const pmlEta & Eta,const std::string Sign):basisMatFunc("JpmlRexp"+Sign),eta(Eta),plus(Sign=="Plus"){if(Sign!="Plus" and Sign !="Minus")ABORT("specify sign as 'Plus' or 'Minus', is: "+Sign);}
public:
    complex<double> operator()(complex<double> z) const {
        if(plus)return exp(complex<double>(0, 2.)*eta(z)*(z-eta.r0upper()))/z;
        else    return exp(complex<double>(0,-2.)*eta(z)*(z-eta.r0upper()))/z;}
    bool operator==(const basisMatFunc & Other) const {
        if(name!=Other.name)return false;
        const pmlRexpSig * o=dynamic_cast<const pmlRexpSig*>(&Other);
        if(o==0)ABORT("two different basisMatFunc instances with same name: "+o->name);
        return o->eta==eta and o->plus==plus;}
};
class pmlRexpRoot:public basisMatFunc{
    friend class ComplexScaling;
    pmlEta eta;
    bool plus;
    pmlRexpRoot(const pmlEta & Eta,const std::string Sign ):basisMatFunc("JpmlRoot"+Sign),eta(Eta),plus(Sign=="Plus"){if(Sign!="Plus" and Sign !="Minus")ABORT("specify sign as 'Plus' or 'Minus', is: "+Sign);}
public:
    complex<double> operator()(complex<double> z) const {
        if(plus)return exp(complex<double>(0, 1.)*eta(z)*(sqrt(z)-sqrt(eta.r0upper())))/z;
        else    return exp(complex<double>(0,-1.)*eta(z)*(sqrt(z)-sqrt(eta.r0upper())))/z;}
    bool operator==(const basisMatFunc & Other) const {
        if(name!=Other.name)return false;
        const pmlRexpRoot * o=dynamic_cast<const pmlRexpRoot*>(&Other);
        if(o==0)ABORT("two different basisMatFunc instances with same name: "+o->name);
        return o->eta==eta and o->plus==plus;}
};

/// @endcond

/// pick complex scaling parameters for AxisName from input
/// if AxisName not found in input list, default to ECS with infinite unscaled region
ComplexScaling::ComplexScaling(ReadInput & in, string AxisName)
    :eta(1.),_r0upper(DBL_MAX),_r0lower(-DBL_MAX),kind("ECS")
{
    double theta;
    bool unitary;
    double rGauge=0;
    if(Algebra::isAlgebra("Rg"))rGauge=real(Algebra("Rg").val(0.));
    in.texdocuCategoryAdd("Absorption","axis,theta,upper,lower,kind","","05,09,23");
    for(unsigned int l=1;;l++){
        in.read("Absorption","axis",axis,"NONE","absorption parameters for this axis",l)
                .texdocu(R"tex(
                         must macth one of the infinite, unbound axes in \nameref{docu:Axis:name}
                         )tex");
        if(axis=="NONE" and not in.noInputFile())break;
        ComplexScaling::names.push_back(axis);
        if(axis!=AxisName and not in.noInputFile())continue;

        in.read("Absorption","kind",kind,"ECS","complex scaling kind",l,"","PML")
                .texdocu(R"tex(
                         Defaults to ECS (exterior complex scaling), this is the only option that is maintained at present.
                         )tex");
        in.read("Absorption","theta",theta,"0.1","complex parameter eta=exp(i theta)",l,"","[0,1.571]")
                .texdocu(R"tex(
                         Only in ECS, allowed values $\in [0,\pi/2]$, good values in ranges $0.1\sim0.5$, needs trying for new systems.
                         Larger values can be more efficient. In case of instability, try smaller values.
                         )tex");
        in.read("Absorption","upper",_r0upper,tools::str(DBL_MAX),"upper end of inner region",l)
                .texdocu(R"tex(
                         Main complex scaling radius $R_0$ for the axis (for any coordinate, e.g. \lcode{X,Rn,Rho})
                         )tex");
        in.read("Absorption","lower",_r0lower,tools::str(-_r0upper),"lower end of inner region",l)
                .texdocu(R"tex(
                         Negative complex scaling radius for coordinates $(-\infty,\infty)$: \lcode{X,Y,Z}.
                         Defaults to -\nameref{docu:Absorption:upper}.
                         )tex");
        in.read("Absorption","unitary",unitary,"false","interprete as the equivalent unitary transformation",l)
                .texdocu(R"tex([DEVELOPER])tex");

        if(rGauge>0. and _r0upper<rGauge and rGauge<DBL_MAX/2)ABORT("choose radial complex scaling radius >= gauge switching radius");

        if(in.noInputFile())return;

        if(kind=="ECS" or kind=="modECS"){
            if(unitary)eta=exp(theta);
            else eta=complex<double>(cos(theta),sin(theta));
        }
        else if(kind=="modPML"){
            if(unitary){
                ABORT("do not know how to do unitary modPML");
            } else {
                eta=1.0+complex<double>(0.0,theta);cout<<"Using eta = 1 + i*"<<eta.imag()<<endl;
            };
        }
        else if(kind=="PML"){
            if(unitary)eta=theta;
            else eta=complex<double>(0.0,theta);
            basisMatFuncSet(new pmlEta(*this));
            basisMatFuncSet(new pmlRexpSig(pmlEta(*this),"Plus"));
            basisMatFuncSet(new pmlRexpSig(pmlEta(*this),"Minus"));
            basisMatFuncSet(new pmlEtaRoot(pmlEta(*this)));
            basisMatFuncSet(new pmlRexpRoot(pmlEta(*this),"Plus"));
            basisMatFuncSet(new pmlRexpRoot(pmlEta(*this),"Minus"));
        }
        else ABORT("Unknown absorbing boundary method '"+kind+"'");
    }
    axis=AxisName;
}

std::string ComplexScaling::strDefinition() const {
    //    std::string axis;
    //    std::complex<double> eta;      ///< complex scaling angle
    //    double _r0upper,_r0lower;        ///< upper/lower scaling radius
    //    std::string kind;              ///< kinds: ECS, PML (may be extended)
    Str s("ComSca: "+kind,",");
    s=s+axis+eta.real()+eta.imag()+_r0lower+_r0upper;
    return s;
}

ComplexScaling::ComplexScaling(string Def){
    if(Def.find("ComSca:")!=0)DEVABORT("not a ComSca definition string "+Def);
    std::vector<std::string> par(tools::splitString(Def.substr(8),','));
    kind=tools::cropString(par[0]);
    axis=tools::cropString(par[1]);
    eta=std::complex<double>(tools::string_to_double(par[2]),tools::string_to_double(par[3]));
    _r0lower=tools::string_to_double(par[4]);
    _r0upper=tools::string_to_double(par[5]);
}

// return local eta
complex<double> ComplexScaling::etaX(double X) const {
    if(kind=="ECS" or kind=="PML"){
        if  (tools::doubleBelow(X,_r0lower,_r0upper))
            return eta;
        else if(tools::doubleAbove(X,_r0lower,_r0upper))
            return eta;
        else
            if     (kind=="ECS")return 1.; // unscaled element
            else if(kind=="PML")return 0.; // unscaled element
            else
                ABORT("etaX not fully defined for ComplexScaling.kind="+kind);
    }
    else if(kind=="modECS"){
        if     (X > _r0upper)return      eta;
        else if(X < _r0lower)return conj(eta);
        else return 1.0; // unscaled element
    }
    else if(kind=="PML"){
        if     (X > _r0upper or X < _r0lower)return eta;
        else return 0.0; // unscaled element
    }
    else if(kind=="UNSCALED")return 1.;
    else {
        ABORT("etaX not defined for comsca.kind="+kind);
        return 0.;
    }
}

ComplexScaling::ComplexScaling(std::string Axis, double theta, double r0upper, double r0lower, std::string Kind)
    :axis(Axis),_r0upper(r0upper),_r0lower(r0lower),kind(Kind) {eta=complex<double>(cos(theta),sin(theta));}

void ComplexScaling::coordinates(UseMatrix &z) const{
    for( unsigned int k=0;k<z.size();k++){
        complex<double> x=z(k).complex();
        if(imag(x)!=0.)ABORT("cannot complex scale, coordinate is complex: "+tools::str(x));
        z(k)=xScaled(x.real());
    }
}
bool ComplexScaling::operator==(const ComplexScaling & B) const {
    if(eta==1.)return eta==B.eta;
    return (axis==B.axis)&&(eta==B.eta)&&(kind==B.kind)&&(_r0lower==B._r0lower)&&(_r0upper==B._r0upper);
}

complex<double> ComplexScaling::xScaled(double X) const{
    if(kind=="ECS" or kind=="modPML"){
        if(X>_r0upper)return _r0upper+eta*(X-_r0upper);
        if(X<_r0lower)return _r0lower+eta*(X-_r0lower);
        else return complex<double>(X);
    }
    else if(kind=="PML" or kind=="UNSCALED")
        return complex<double>(X);
    else
        ABORT("coordinate transformation not defined for ComplexScaling.kind="+kind);
}

string ComplexScaling::str() const{
    if(eta==1.)return "";
    return kind+":("+tools::str(eta.real(),3)+","+tools::str(eta.imag(),3)+")["+tools::str(_r0lower,3)+","+tools::str(_r0upper,3)+"]";
}
