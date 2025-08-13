// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1 // needed, as std::sph_bessel not implemented in c++11 by default

#include "algebra.h"
#include <cmath>
#include "constants.h"
#include "printOutput.h"
#ifdef _USE_BOOST_
#include <boost/math/special_functions.hpp>
#endif
#include "integrate.h"
#include "readInput.h"
#include "asciiFile.h"
#include "algebraMulti.h"

using namespace std;

const std::string Algebra::Infty="infty";
string Algebra::failures;
std::vector<const Algebra*> Algebra::integrand; // holds Algebra's for nested integration, .back() is latest

std::map<std::string,std::complex<double> >Algebra::specialConstants;
std::map<std::string,std::complex<double> >Algebra::updatableConstants;
std::map<std::string,Algebra::algebraFactory> Algebra::externalFactory;

bool Algebra::isSpecialConstant(std::string Term){
    specialConstants["i"]=complex<double>(0.,1.);
    specialConstants["pi"]=math::pi;
    specialConstants["hbar"]=physics::h_bar;
    specialConstants["alfaFine"]=physics::a_finestructure;
    specialConstants["infty"]=DBL_MAX;

    return specialConstants.count(Term);
}

bool Algebra::isUpdatableConstant(std::string Term){return updatableConstants.count(Term);}


void Algebra::addExternalFactory(algebraFactory Factory,std::string Name){
    if(Name=="from FactoryName")
        Name=std::shared_ptr<const Algebra>(Factory(""))->definition();
    externalFactory[Name]=Factory;
}

void Algebra::addSpecialConstant(string Name, std::complex<double> Value){
    if(isAlgebra(Name)){
        if(isSpecialConstant(Name) and specialConstants[Name]==Value)return;
        ABORT( "\nillegal constant name: "+Name+" with value="+tools::str(Value,14)
               +"\nconflicts with function name "
               +"\nor with the previously defined constants: "+listConstants());
    }
    specialConstants[Name]=Value;
    specialConstants["-"+Name]=-Value;
}

void Algebra::addUpdatableConstant(string Name, std::complex<double> Value){
    updatableConstants[Name]=Value;
    updatableConstants["-"+Name]=-Value;
}

void Algebra::readConstants(ReadInput &Inp){
    string name;
    double value;

    for(unsigned int line=1;;line++){
        Inp.read("Constant","name",name,"BLANK","name of constant for use in algebra",line)
                .texdocu(R"tex(
                         defines name to be used in \nameref{sec:class:Algebra} type expressions
                         )tex");
        if(name=="BLANK")break;
        Inp.read("Constant","value",value,ReadInput::noDefault,"value of constant corresponding to name",line)
                .texdocu(R"tex(
                         value associated with \nameref{docu:Constant:name}
                         )tex");
        addSpecialConstant(name,value);
        listConstants();
    }
}


bool Algebra::isConstant(std::string Term)
{
    if(Term.find("const[")==0)return true;
    return  Term.substr(0,Term.find("[")).find_first_not_of(" +-.0123456789e,")>Term.substr(0,Term.find("[")).length();
}

static std::map<std::string,std::string> standardAlgebras;
std::string Algebra::listStandard(){
    Algebra("Q");
    std::string res;
    for(auto s: standardAlgebras){
        res+=s.first+std::string(15-s.first.length(),' ')+s.second+"\n";
    }
#ifdef _DEVELOP_
    res+="!!! DEVELOPER: check Algebra::factory for more !!!\n";
#endif
    return res;
}
const Algebra* Algebra::factory(const string Term, bool Add){
    const Algebra* amult=AlgebraMulti::factory(Term);
    if(amult)return amult;

    if(isConstant(Term)){return new AlgebraConstant(Term);}
    if(isSpecialConstant(Term))return new AlgebraConstant(Term,specialConstants[Term]);
    if(isUpdatableConstant(Term))return new AlgebraConstantUpdatable(Term);

    bool isStandard=Term=="Q" or Term.find("valUp")==0 or Term.find("derUp")==0;
    // this is for info and needs to duplicate the choices below (until more elegant solution is found)
    standardAlgebras={
        {"sqrt","square-root"}
        ,{"const","const[val](Q)...constant value (also specify as 'val' w/o const[..]"}
        ,{"sin","sin(Q)...sine"}
        ,{"cos","cos(Q)...cosine"}
        ,{"exp","exp[al](Q)...exp(al*Q), al a legal float string"}
        ,{"ExpI","ExpI[al](Q)...exp(i*al*Q), al a legal float string "}
        ,{"pow","pow[n](Q)...(integer) power Q^n"}
        ,{"chi","chi[a,b](Q)...characteristic function on [a,b], a,b legal float strings"}
        ,{"trunc","trunc[a,b]...=1 on [-a,a], =0 outside [-b,b], differentiably smooth between (3rd order polynomial)"}
        ,{"spherBessel","spherBessel[l](Q)...spherical Bessel function for l"}
        ,{"Morse","Morse[stiff,minpos](Q)...Morse potential (see AlgebraMorse for details)"}
    };
    for(auto s: standardAlgebras){
        if(isStandard)break;
        isStandard=Term.find(s.first)==0;
    }
    if(isStandard){
        // possible solution: introduce Algebra* AlgebraSqrt::try(Term) which returns 0 if fail
        if(Term=="Q")return new AlgebraQ();
        if(Term.substr(0,5)=="const[")return new AlgebraConstant(Term);
        if(Term.substr(0,5)=="sqrt(")return new AlgebraSqrt(Term);
        if(Term.substr(0,4)=="sin(")return new AlgebraSin(Term);
        if(Term.substr(0,4)=="cos(")return new AlgebraCos(Term);
        if(Term.substr(0,4)=="exp(")return new AlgebraExp(Term);
        if(Term.substr(0,4)=="chi[")return new AlgebraChi(Term);
        if(Term.substr(0,6)=="trunc[")return new AlgebraTrunc(Term);
        if(Term.substr(0,7)=="valUp1[")return new AlgebraValUp1(Term);
        if(Term.substr(0,7)=="derUp1[")return new AlgebraDerUp1(Term);
        if(Term.substr(0,4)=="pow[")return new AlgebraPow(Term);
        if(Term.substr(0,12)=="spherBessel[")return new AlgebraSpherBessel(Term);
        if(Term.substr(0,5)=="ExpI[")return new AlgebraExpI(Term);
        if(Term.substr(0,6)=="Morse[")return new AlgebraMorse(Term);
    }
    string name=Term.substr(0,Term.find_first_of("(["));
    if(externalFactory.count(name)==1)return externalFactory[name](Term);

    const Algebra* res(0);
    // last resort: try all external factories (if any)
    for(auto f: externalFactory){
        if(f.first!="" and (res=f.second(Term)))return res;
    }

    addFailure("UNDEFINED: "+Term);
    return 0;
}

Algebra::~Algebra(){
    for(unsigned int k=0;k<A.size();k++)delete A[k];
    for(unsigned int k=0;k<I.size();k++)delete I[k];
}

bool Algebra::isAlgebra(string Definition){
    Algebra a(Definition);
    return a.isAlgebra();
}

bool Algebra::isAlgebra() const{
    for(unsigned int k=0;k<A.size();k++)
        if(A[k]==0 or not A[k]->isAlgebra())return false;
    for(unsigned int k=0;k<I.size();k++)
        if(I[k]==0 or not I[k]->isAlgebra())return false;
    return true;
}

std::string Algebra::str() const {
    if(not isAlgebra())return failures;
    std::string op=add?"+":"*";
    std::string res;
    res+=Str("","")+A.size()+"/"+I.size();
    for(const auto a: A){
        res+=op+": "+a->_definition+"\n";
        res+=a->str();
    }
    op=op=="+"?"-":"/";
    for(const auto a: I){
        res+=op+": "+a->_definition+"\n";
        res+=a->str();
    }
    return res;
}

bool Algebra::isAlgebraOfConsts() const{
    for(unsigned int k=0;k<A.size();k++)
        if(A[k]==0 or not A[k]->isAlgebraOfConsts())return false;
    for(unsigned int k=0;k<I.size();k++)
        if(I[k]==0 or not I[k]->isAlgebraOfConsts())return false;
    return (A.size()+I.size())>0 or dynamic_cast<const AlgebraConstant *>(this)!=0;
}

bool Algebra::isAlgebraOfUpdatables() const{
    bool hasUpd=false;
    for(unsigned int k=0;k<A.size();k++)
        if(A[k] and A[k]->isAlgebraOfUpdatables())hasUpd=true;
    for(unsigned int k=0;k<I.size();k++)
        if(I[k] and I[k]->isAlgebraOfUpdatables())hasUpd=true;
    if((A.size()+I.size())==0)hasUpd=dynamic_cast<const AlgebraConstantUpdatable*>(this);
    return hasUpd;
}

double Algebra::signChangeQ(double Low, double Up, double ValueLow, double ValueUp, double Epsilon) const {
    if(signbit(ValueLow)!=signbit(ValueUp))DEVABORT("no sign change in interval");
    double center((Low+Up)/2);
    if(Low-Up<=Epsilon)return center;
    double vCenter=val(center).real();
    return signbit(vCenter)!=signbit(ValueLow)?
                signChangeQ(Low,center,ValueLow,vCenter,Epsilon)
              : signChangeQ(center,Up, vCenter,ValueUp, Epsilon);
}

std::vector<double> Algebra::zeros(double Low, double Up, double Delta, double Epsilon) const {
    if(Up-Low<=Delta){
        double vLow=val(Low).real();
        double vUp=val(Up).real();
        if(signbit(vLow)!=signbit(vUp))
            std::vector<double>(signChangeQ(Low,Up,vLow,vUp,Epsilon));
    }
    std::vector<double>zero(zeros(Low,(Low+Up)/2,Delta));
    std::vector<double>zer1(zeros((Low+Up)/2,Up,Delta));
    zero.insert(zero.end(),zer1.begin(),zer1.end());
    return zero;
}

std::vector<std::complex<double>> Algebra::nonAnalyticQ() const {
    std::vector<std::complex<double>> singQ;
    for(auto a: A){
        auto sing(a->nonAnalyticQ());
        singQ.insert(singQ.end(),sing.begin(),sing.end());
    }
    for(auto a: I){
        auto sing(a->nonAnalyticQ());
        singQ.insert(singQ.end(),sing.begin(),sing.end());
    }
    return singQ;
}

std::vector<double> Algebra::poles(double Low, double Up, double Delta, double Epsilon) const {
    std::vector<double> singQ;
    for(auto a: A){
        auto sing(a->poles(Low,Up,Delta,Epsilon));
        singQ.insert(singQ.end(),sing.begin(),sing.end());
    }
    for(auto a: I){
        std::vector<double> sing;
        if(add)sing=a->poles(Low,Up,Delta,Epsilon);
        else   sing=a->zeros(Low,Up,Delta,Epsilon);
        singQ.insert(singQ.end(),sing.begin(),sing.end());
    }
    return singQ;
}

std::vector<std::complex<double>> AlgebraTrunc::nonAnalyticQ() const{
    std::vector<std::complex<double>> sing({-from,-to,from,to});
    return sing;
}

complex<double> Algebra::val(const std::complex<double> Q) const {
    complex<double> res,ires;
    if(add){
        res=0.;
        ires=0;
        for(unsigned int k=0;k<A.size();k++) res+=A[k]->val(Q);
        for(unsigned int k=0;k<I.size();k++)ires+=I[k]->val(Q);
        return res-ires;
    }else{
        res=1.;
        ires=1.;

        for(unsigned int k=0;k<A.size();k++){
            complex<double> fac=A[k]->val(Q);
            if(fac==0.){res=0.;break;}
            res*=fac;
        }
        for(unsigned int k=0;k<I.size();k++){
            complex<double> fac=I[k]->val(Q);
            if(fac==0.)ABORT(Str("denominator zero at Q=")+Q+"for"+_definition+" from factor "+I[k]->_definition);
            ires*=fac;
        }
        return res/ires;
    }
}

static complex<double> Func(const vector<double>&X)
{
    if(X.size()!=1)ABORT("need vector length =1");
    return Algebra::integrand.back()->val(complex<double>(X[0]));
}
complex<double> Algebra::integral(const std::complex<double> Q0, const std::complex<double> Q1) const{
    if(Q0==Q1)return 0.;

    integrand.push_back(this);
    Int integ;

    if(Q0.imag()!=0 or Q1.imag()!=0)ABORT("integration over complex contours not implemented");
    vector<vector<double> > vol(vector<vector<double> >(1));
    std::complex<double> res;
    if(Q0.real()<Q1.real()){
        vol[0].push_back(Q0.real());
        vol[0].push_back(Q1.real());
        res= integ.recursive(vol,Func);
    } else {
        vol[0].push_back(Q1.real());
        vol[0].push_back(Q0.real());
        res=-integ.recursive(vol,Func);
    }
    integrand.pop_back();
    return res;
}


bool allowedCharacters(string Def){
    for(size_t k=0;k<Def.length();k++){
        if(not isalnum(Def[k]) and string("+-*/[](),. ").find(Def[k])==string::npos)return false;
    }
    return true;
}

Algebra::Algebra(const string Definition, const bool Add):add(Add),_definition(Definition){

    if(failures=="")failures=_definition;

    for(size_t k=0;k<Definition.length();k++){
        if(not isalnum(Definition[k]) and string("+-*/[](),. ").find(Definition[k])==string::npos){
            addFailure("ILLEGAL CHARACTER '"+Definition.substr(k,1)+"' "+_definition);return;
        }
    }
    if(tools::subStringCount(Definition,"(")!=tools::subStringCount(Definition,")")){
        addFailure("UNBALANCED PARENTHESIS "+_definition); return;
    }
    if(tools::findFirstOutsideBrackets(Definition,",","([",")]")!=string::npos){
        addFailure("COMMA "+_definition); return;
    }

    string def1=tools::cropString(_definition);

    // split into terms
    vector<string> terms,kind;
    if(add)tools::splitString(def1,"+-",terms,kind,"([",")]");
    else   tools::splitString(def1,"*/",terms,kind,"([",")]");

    // undo split at exponents: ends in 'e' preceded by digit or '.'
    if(add){
        for(size_t k=1;k<terms.size();k++){
            string s=terms[k-1];
            if(s.length()>1 and s[s.length()-1]=='e' and string("01234567890.").find(s[s.length()-2])!=string::npos){
                terms[k-1]+=kind[k]+terms[k];
                kind.erase(kind.begin()+k);
                terms.erase(terms.begin()+k);
            }
        }
    }

    // multiple terms
    if(terms.size()>1 or add){
        // collect direct and inverse terms
        for(unsigned int k=0;k<terms.size();k++){
            if(terms[k]==""){addFailure("EMPTY TERM AT "+_definition);return;}
            if(kind[k][0]=='-' or kind[k][0]=='/')
                I.push_back(new Algebra(terms[k],not add));
            else
                A.push_back(new Algebra(terms[k],not add));
        }
    }

    // single term
    else{
        string def0=tools::cropString(terms[0]);
        if(def0[0]!='('){
            // no further brackets at this term
            if(kind[0][0]=='/')I.push_back(factory(def0,true));
            else               A.push_back(factory(def0,true));
        }
        else {
            // remove outermost pair of (...)
            if(def0[def0.length()-1]==')')def0=def0.substr(1,def0.length()-2);
            A.push_back(new Algebra(def0,true));
        }
    }
}

string Algebra::argument(string Term) const {
    string s=tools::cropString(Term);
    if(s.find(_definition+"(")!=0 or s[s.length()-1]!=')'){
        if(tools::findOutsideBrackets(true,_definition,"(","[","]")==string::npos
                and tools::findOutsideBrackets(true,_definition,")","[","]")==string::npos)return "Q"; // assume default argument
        ABORT("function \""+_definition+"(...)\" must have argument, found "+Term);
    }
    return s.substr(_definition.length()+1,s.length()-_definition.length()-2);
}

std::complex<double> Algebra::constantValue(string Term){
    Algebra a(Term);
    if(not a.isAlgebraOfConsts())ABORT("not a constant expression: "+Term);
    return a.val(0.);
}

double Algebra::realConstant(std::string Term){
    std::complex<double> res=constantValue(Term);
    if(res.imag()!=0.)ABORT("not a real number: "+Term);
    return res.real();
}

int Algebra::integerConstant(std::string Term){
    double res=Algebra::realConstant(Term);
    if(res-int(res)!=0.)ABORT("not an integer number: "+Term);
    return int(res);
}

double Algebra::smooth3rdDegree(double x, double a, double b, int der){
    double fac=-1./(b-a);
    double q=(b-x)/(b-a);
    switch(der){
    case(0): return q*q*(3.-2.*q);
    case(1): return q*(6.-6.*q)*fac;
    case(2): return (6.-12.*q)*fac*fac;
    default: ABORT("n'th derivative of trunc not defined, n="+tools::str(der));
    }
    return 0.;
}

double Algebra::valUp1(double x, double a, double b){
    double q=(x-a)/(b-a);
    return q*(2.-q);
}

double Algebra::derUp1(double x, double a, double b){
    double q=(x-a)/(b-a);
    return q*(q-1)*(b-a);
}

double Algebra::getParameter(unsigned int k, string Definition){
    //    if(Definition.find("]")+1!=Definition.length())ABORT("for parameters, must have ] at the end of the definition, is: "+Definition);
    std::vector<std::string> pars=tools::splitString(tools::stringInBetween(Definition,"[","]"),',');
    if(pars.size()<k+1)ABORT("need at least "+tools::str(k+1)+" parameters, have: "+Definition);
    //    if(isConstant(pars[k]))return std::real(AlgebraConstant(pars[k]).val(0.));
    if(isSpecialConstant(pars[k]))return std::real(specialConstants[pars[k]]);
    Algebra a(pars[k]);
    if(a.isAlgebraOfConsts())return std::real(a.val(0.));
    if(a.isAlgebra()){
        if(Definition.find("Q")!=string::npos)
            ABORT("a non-constant function parameter seems to be used: "+Definition+" (must not contain variable Q)");
        return std::real(a.val(0.));
    }
    else
        ABORT("not an algebraic expression: "+Definition+", syntax errors: "+a.failures);
    ABORT("non-constant function parameter in "+Definition+"\navailable constants: \n"+listConstants());
    return 0.;
}

complex<double> AlgebraSpherBessel::val(std::complex<double> k) const {
    switch (der) {
    case 0:// value
        // assume real k for sph_bessel; does not work with complex<double>
#ifdef _USE_BOOST_
        return pow(2./M_PI,.5)*pow(std::complex<double>(0,-1.),order)*
                boost::math::sph_bessel(order, k.real()*surface);
#else
        return pow(2./M_PI,.5)*pow(std::complex<double>(0,-1.),order)*
                std::sph_bessel(order, k.real()*surface);
#endif

    case 1: // derivative; use recurrence relation
#ifdef _USE_BOOST_
        return  pow(2./M_PI,.5)*pow(std::complex<double>(0,-1.),order)
                *(boost::math::sph_bessel(order, k.real()*surface)/surface*order
                  -boost::math::sph_bessel(order+1, k.real()*surface)*k.real());
#else
        return  pow(2./M_PI,.5)*pow(std::complex<double>(0,-1.),order)
                *(std::sph_bessel(order, k.real()*surface)/surface*order
                  -std::sph_bessel(order+1, k.real()*surface)*k.real());
#endif
    default: ABORT("der can only be 0 or 1, is: "+tools::str(der));
    }
    return 0.;
}

void Algebra::Test(){
    vector<complex<double> > result;
    vector<string> def;
    def.push_back("-0.5");     result.push_back(-0.5);
    def.push_back("2*Q-sqrt(3.*Q*Q+Q/(1./Q))");     result.push_back(0.);
    def.push_back("sin(Q)*sin(Q)+cos(Q)*cos(Q)-1");     result.push_back(0.);
    def.push_back("exp(Q)/exp(-Q)-exp(2*Q)");     result.push_back(0.);
    def.push_back("chi[0.,1.](0.3*Q)");     result.push_back(1.);
    def.push_back("chi[-1.,1.](-Q)");     result.push_back(0.);
    def.push_back("trunc[3.,7.](2.)");     result.push_back(1.);
    def.push_back("trunc[3.,7.](7.)");     result.push_back(0.);
    def.push_back("trunc[3.,7.](6.)");     result.push_back(0.25*0.25*(3.-0.25*2));
    def.push_back("pow[1.5](pow[2](Q))-Q*Q*Q"); result.push_back(0.);

    // deliberatly malformed strings
    def.push_back("Exp(Q)/exp(-Q)-exp(2*Q)");result.push_back(0.);
    def.push_back("(exp(Q)/exp(-Q)-exp(2*Q)");result.push_back(0.);
    def.push_back("(exp(Q)/exp(-Q)-exp(2*Q))");result.push_back(0.);
    def.push_back("(exp(Q)*/exp(-Q)-exp(2*Q))");result.push_back(0.);
    def.push_back("(exp(Q)//exp(-Q)-exp(2*Q))");result.push_back(0.);

    vector<const Algebra*> a;
    for(unsigned int k=0;k<def.size();k++){
        a.push_back(new Algebra(def[k]));
        if(not a.back()->isAlgebra()){
            cout<<failures<<endl;
            a.pop_back();
        }
        failures="";
    }

    complex<double>arg=1.414;
    for(unsigned int k=0;k<a.size();k++){
        if(abs(a[k]->val(arg)-result[k])>1.e-12)
            cout<<" OK: "+a[k]->_definition+" = "<<result[k]<<" at Q="<<arg<<endl;
        else
            ABORT("failed to evaluate "+a[k]->_definition+" at Q="+tools::str(arg)+"=?="+tools::str(a[k]->val(arg)));

    }
    cout<<"Algebra::test() passed"<<endl;

    ofstream plot;
    plot.open("algebraPlots");
    vector<Algebra*> plt;
    plt.push_back(new Algebra("trunc[3.,7.](Q)"));
    plt.push_back(new Algebra("pow[0.3](pow[2](Q))"));
    plt.push_back(new Algebra("chi[0.,1.](0.3*Q)"));
    plot<<"#";
    cout<<"*** plots \"algebraPlots\":";
    for(unsigned int k=0;k<plt.size();k++){
        plot<<",  "<<plt[k]->_definition;
        cout<<" "<<plt[k]->_definition+",";
    }
    cout<<endl;
    plot<<endl;
    for (double x=-10.;x<=10.;x+=0.1){
        plot<<x;
        for(unsigned int k=0;k<plt.size();k++)plot<<", "<<plt[k]->val(x).real();
        plot<<endl;
    }
}

void Algebra::plot(string File, double From, double To, int Points) const{
    if(not isAlgebra())ABORT("malformed algebra string: "+_definition+" "+failures);

    vector<vector<double> >cols(2);
    for(int k=0;k<Points;k++){
        cols[0].push_back(From+k*(To-From)/(Points-1));
        cols[1].push_back(val(cols[0][k]).real());
    }

    vector<string> coms(1,_definition);
    AsciiFile f(File);
    f.writeComments(coms);
    f.writeCols(cols);
    PrintOutput::message(Str("Algebra")+_definition+"plot on"+File);
}

std::complex<double> AlgebraExpI::val(const std::complex<double> Q) const
{
    return pow(std::complex<double>(0., 1)*Q, der)*exp(std::complex<double>(0., surface)*Q)/sqrt(2.*M_PI);
}
