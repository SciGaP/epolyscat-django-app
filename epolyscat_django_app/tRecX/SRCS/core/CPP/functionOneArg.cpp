// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "functionOneArg.h"
#include "constants.h"
#include "userFunctions.h"
#include "orthopol.h"
#include "pulse.h"
//#include "basisMat.h"

using namespace std;

/** @name Standard function classes
 * constructors private, only friend class FunctionOneArg  can use these
*/

/// get pointer from list or try to add to list if missing
FunctionOneArg* FunctionOneArg::get(std::string Name, bool Abort){

    if(not tools::hasKey(list,Name)){

        // get all available functions
        standard();

        // check for duplicates
        for(unsigned int k=0;k<avail.size();k++)
            for(unsigned int l=0;l<k;l++)
                if(avail[k]->name==avail[l]->name)
                    ABORT("\n duplicate function name: "+avail[k]->name+"\n check FunctionOneArg derived classes");

        // select suitable
        string fName=Name.substr(0,Name.find("["));
        for(unsigned int k=0;k<avail.size();k++){
            if(fName==avail[k]->name)list[Name]=avail[k];
        }

        if(Abort and not tools::hasKey(list,Name)){
                string availNames;
                for(unsigned int k=0;k<avail.size();k++)availNames+="\n "+avail[k]->name;
                ABORT("\n function '"+fName+"' ('"+Name+"') not defined\n Available names:"+availNames+"\n");
        }

        // remove unneeded
        for(unsigned int k=0;k<avail.size();k++)
            if(avail[k]->name!=fName)delete avail[k];
        avail.clear();

        // no such function, but not aborted:
        if(not tools::hasKey(list,Name))return 0;

        // set parameters
        vector<string> sPars = tools::splitString(tools::stringInBetween(Name,"[","]"),',');
        if(sPars.size()==1 and sPars[0]=="t")sPars.clear(); // paramter "t" will be considered as a (time) argument
        vector<double> fPars;
        for(unsigned int i=0;i<sPars.size();i++)fPars.push_back(tools::string_to_double(sPars[i]));
        list[Name]->pars=fPars;
        list[Name]->defaults(); // supplement with defaults (if needed)
    }
    return list[Name];
}

/// @cond DEV
///@{

/// pulse[i]... pulse as defined on input, i=0,1,2: z,x,y-component
class FunctionOneArg::FunctionPulse:public FunctionOneArg{
    friend class FunctionOneArg;
    Pulse p;
    FunctionPulse():FunctionOneArg("pulse"){}
    void defaults(){
        p=Pulse::current;
        if(p.sing.size()==0)ABORT("no pulse defined, call Pulse::read(Inp,\"Pulse\")");
        p.output("PULSE");
        if(pars.size()!=1 or pars[0]>2)ABORT("must define pulse component 0,1,2");
    }
    std::complex<double> val(double Q) const{
        return p.Apot(Q,(unsigned int)pars[0]);}
};

/// const[p0]... constant with value p0
class FunctionOneArg::FunctionConst: public FunctionOneArg{
    friend class FunctionOneArg;
    FunctionConst():FunctionOneArg("const"){}
    void defaults() {if(pars.size()!=1)ABORT(name+"needs exactly one parameter");}
    complex<double> val(double Q) const {return pars[0];}
    complex<double> val(complex<double> Q) const {return pars[0];}
};

/// delta[p0]...\f$ \delta(q-p_0)\f$
class FunctionOneArg::FunctionDelta: public FunctionOneArg{
    friend class FunctionOneArg;
    FunctionDelta():FunctionOneArg("delta"){}
    void defaults() {if(pars.size()!=1)ABORT(name+"needs exactly one parameter");}
    complex<double> val(double Q) const {if(Q==pars[0]) return 1; return 0.;} ///< returns \f$\cos( par[0] q )\f$
    complex<double> val(complex<double> Q) const {if(Q==pars[0]) return 1; return 0.;} ///< returns \f$\cos( par[0] q )\f$
};

/// cos[p0]...\f$ \cos(p_0\times q)\f$
class FunctionOneArg::FunctionCos: public FunctionOneArg{
    friend class FunctionOneArg;
    FunctionCos():FunctionOneArg("cos"){}
    void defaults() {if(pars.size()!=1)ABORT(name+"needs exactly one parameter");}
    complex<double> val(double Q) const {return cos(pars[0]*Q);} ///< returns \f$\cos( par[0] q )\f$
};

/// sin[p0]...\f$ \sin(p_0\times q)\f$
class FunctionOneArg::FunctionSin:public FunctionOneArg{
    friend class FunctionOneArg;
    FunctionSin():FunctionOneArg("sin"){}
    void defaults() {if(pars.size()!=1)ABORT(name+"needs exactly one parameter");}
    complex<double> val(double Q) const {return sin(pars[0]*Q);} ///< returns \f$\cos( par[0] q )\f$
};

/// gauss[p0,p1,p2]...\f$ \exp[- (q-p_1)^2/p_0^2]\f$ truncated at \f$ |q|>p_2 \f$ (defaults: p0=1,p1=0,p2=infty)
class FunctionOneArg::FunctionGauss:public FunctionOneArg{
    friend class FunctionOneArg;
    FunctionGauss():FunctionOneArg("gauss"){}
    void defaults() {
        if(pars.size()==0)ABORT(name+"needs at least one parameter");
        pars.resize(3,DBL_MAX);
        if(pars[0]==DBL_MAX)pars[0]=1.;
        if(pars[1]==DBL_MAX)pars[1]=0.;
    }
    complex<double> val(double Q) const
    {   if(abs(Q)>pars[2])return 0;
        return exp(-pow((Q-pars[1])/pars[0],2));
    }
};

/// rGauss[p0,p1,p2]...\f$ q\exp[- (q-p_1)^2/p_0^2]\f$ truncated at \f$ |q|>p_2 \f$ (defaults: p0=1,p1=0,p2=infty)
class FunctionOneArg::FunctionRGauss:public FunctionOneArg{
    friend class FunctionOneArg;
    FunctionRGauss():FunctionOneArg("rGauss"){}
    void defaults() {
        pars.resize(3,DBL_MAX);
        if(pars[0]==DBL_MAX)pars[0]=1.;
        if(pars[1]==DBL_MAX)pars[1]=0.;
    }
    complex<double> val(double Q) const
    {   if(abs(Q)>pars[2])return 0;
        return Q*exp(-pow((Q-pars[1])/pars[0],2));
    }
};

/// cos2[p0,p1,p2]...\f$ \cos^2[(q-p_1)/p_0] \f$ truncated at \f$ |q-p_1|>p_2 \f$ (defaults: \f$ p_0=1,p_1=0,p_2=\pi p_0/2\f$)
class FunctionOneArg::FunctionCos2:public FunctionOneArg{
    friend class FunctionOneArg;
    FunctionCos2():FunctionOneArg("cos2"){}
    void defaults() {
        pars.resize(3,DBL_MAX);
        if(pars[0]==DBL_MAX)pars[0]=1.;
        if(pars[1]==DBL_MAX)pars[1]=0.;
        if(pars[2]==DBL_MAX)pars[2]=pars[0]*math::pi/2;
    }
    complex<double> val(double Q) const
    {   if(abs(Q-pars[1])>pars[2])return 0;
        return pow(cos(Q-pars[1])*pars[0],2);
    }
};

/// legendre[n]...\f$ P_n(q) \f$ Legendre polynomial degree n
class FunctionOneArg::FunctionLegendre:public FunctionOneArg{
    friend class FunctionOneArg;
    FunctionLegendre():FunctionOneArg("legendre"){}
    void defaults() {
        if(pars.size()!=1)ABORT("legendre must have exactly one parameter indicating degree");}
    complex<double> val(double Q) const
    {   OrthogonalPolynomial * p=new OrthogonalLegendre();
        vector<double> v=p->val((int)pars[0]+1,Q);
        return v.back();
    }
};

///// lobatto[n]...\f$ P_n(q) \f$ Lobatto polynomial degree n
//class FunctionOneArg::FunctionLobatto:public FunctionOneArg{
//    friend class FunctionOneArg;
//    FunctionLobatto():FunctionOneArg("lobatto"){}
//    void defaults() {
//        if(pars.size()!=1)ABORT("lobatto must have exactly one parameter indicating degree");}
//    complex<double> val(double Q) const
//    {   OrthogonalPolynomial * p=new OrthogonalLobatto();
//        vector<double> v=p->val((int)pars[0]+1,Q);
//        return v.back();
//    }
//};

///@}

/// @endcond

vector<FunctionOneArg*> FunctionOneArg::avail;
map<string,FunctionOneArg*> FunctionOneArg::list;
void FunctionOneArg::standard(){
    add(new FunctionOneArg::FunctionCos());
    add(new FunctionOneArg::FunctionSin());
    add(new FunctionOneArg::FunctionDelta());
    add(new FunctionOneArg::FunctionGauss());
    add(new FunctionOneArg::FunctionRGauss());
    add(new FunctionOneArg::FunctionPulse());
    add(new FunctionOneArg::FunctionLegendre());
//    add(new FunctionOneArg::FunctionLobatto());
    add(new FunctionOneArg::FunctionConst());
}
