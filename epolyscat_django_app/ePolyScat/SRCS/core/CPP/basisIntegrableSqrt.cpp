#include "basisIntegrableSqrt.h"
#include <complex>

#include "printOutput.h"

#include "basisSetDef.h"
#include "basisIntegrable.h"
#include "basisExpIm.h"
#include "index.h"

int BasisIntegrableSqrt::powerSqrt(std::string Def){
    std::string sSqrt=tools::stringInBetween(Def.substr(0,Def.find("*")),"[","]");
    int _powSqrt(0);
    if(sSqrt=="" or sSqrt=="0")_powSqrt=0;
    else if(sSqrt=="1")_powSqrt=1;
    else ABORT("need power of sqrt 0 or 1, got: "+Def);
    if (_powSqrt<0 or _powSqrt>1)DEVABORT("need power sqrt in [0,1], found: "+Def);
    return _powSqrt;
}

BasisIntegrableSqrt::BasisIntegrableSqrt(const std::string Def):BasisIntegrable(0.,0.){

    PrintOutput::DEVwarning("for serious use, change BasisIntegrableSqrt  by model of BasisIntegrableRpow, or integrate into it",2);

    _basInt=dynamic_cast<const BasisIntegrable*>(BasisAbstract::factory(Def.substr(Def.find("*")+1)));
    if(_basInt==nullptr)ABORT("not an integrable basis after sqrt* :"+Def);
    _lowBound=_basInt->lowBound();
    _upBound=_basInt->upBound();
    if(lowBound()*upBound()<0.)DEVABORT("cannot use sqrt across 0");
    _jac=_basInt->jacobian();
    _comSca=_basInt->complexScaling();
    _powSqrt=powerSqrt(Def);

    _name=Str("sqrt[","")+_powSqrt+"]*"+_basInt->name();
 }

void BasisIntegrableSqrt::valDer(const std::vector<std::complex<double> > &X,
                                 std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const {
    _basInt->valDer(X,Val,Der,ZeroOutside);
    if(_powSqrt==0)return;

    for(size_t k=0,kn=0;k<X.size();k++){
        std::complex<double> sqrt=std::sqrt(std::abs(X[k]));
        std::complex<double> qSqrt2=X[k]==0.?0:0.5/sqrt;
        for(size_t n=0;n<size();n++,kn++){
            Der[kn]=sqrt*Der[kn]+qSqrt2*Val[kn];
            Val[kn]=sqrt*Val[kn];
        }
    }
}

void BasisIntegrableSqrt::resolvePower(BasisSetDef& Def, std::vector<unsigned int> Branch, const std::vector<const Index*> Path){
    if(Def.funcs.substr(0,5)!="sqrt[")return;
    if(Def.funcs.find_first_of("0123456789")==5)return;

    std::string res=Def.funcs;
    std::string powDef=tools::stringInBetween(Def.funcs,"sqrt[","]",true);
    if(powDef.substr(0,4)=="m%2{"){
        int pow=INT_MAX;
        std::string ax=tools::stringInBetween(powDef,"m%2{","}",true);
        for(size_t k=0;k<Path.size();k++){
            if(Path[k]->axisName()==ax and dynamic_cast<const BasisTrigon*>(Path[k]->basis())){
                pow=dynamic_cast<const BasisTrigon*>(Path[k]->basis())->mValueAbs(Branch[k])%2;
                break;
            }
        }
        if(pow==INT_MAX)ABORT("cannot find axis "+ax+" in hierarchy for resolving "+Def.funcs);
        Def.funcs=Str("sqrt[","")+pow+Def.funcs.substr(Def.funcs.find("]"));
    }
    else if(tools::cropString(powDef).find_first_of("0123456789")==std::string::npos)
        ABORT("need format sqrt[m%2{Phi..}]*, do not know how to interprete sqrt power in "+Def.funcs);
}

std::string BasisIntegrableSqrt::inputToStrDefinition(BasisSetDef Def){
    if(Def.funcs.find("sqrt")>2)return "";

    std::string sqrtDef=Def.funcs.substr(0,Def.funcs.find("*"));
    Def.funcs=Def.funcs.substr(Def.funcs.find("*")+1);
    const BasisAbstract* bas=BasisAbstract::factory(Def);
    if(not bas)ABORT("cannot interprete "+Def.str());
    return sqrtDef+"*"+bas->strDefinition();
}
