// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "algebraMulti.h"
//#ifdef _USE_BOOST_
//#include <boost/algorithm/string/replace.hpp>
//#endif
#include "asciiFile.h"

using namespace std;

bool AlgebraMulti::isAlgebraMulti(string Definition){
    std::unique_ptr<const AlgebraMulti> alg(factory(Definition));
    return bool(alg);
}

std::complex<double> AlgebraMulti::val(std::complex<double> Q) const {
    if(_iSingle>_q.size())DEVABORT("single argument use was not constructed");
    _q[_iSingle]=Q;
    return valMulti(_q);
}

const AlgebraMulti* AlgebraMulti::factory(const string Definition){
    if(Definition.find("Sum")==0)     return new AlgebraSum(Definition);
    if(Definition.find("ExternalHO(Eta")==0)return new AlgebraExternalHO(Definition);
    if(Definition.find("CO2Pot")==0)return new AlgebraCO2Pot(Definition);
    if(Definition.find("CO2Nuclear")==0)return new AlgebraCO2Nuclear(Definition);
    if(Definition.find("Numerical[")==0)return new AlgebraNumerical(Definition);
    return 0;
}

AlgebraMulti::AlgebraMulti(std::string Definition, const string Coors):_iSingle(0)
{
    _definition=Definition;
    vector<string> defAlgs=tools::splitString(tools::stringInBetween(Definition,"(",")"),',');
    if(Coors!="*"){
        vector<string> coors=tools::splitString(Coors,',');
        for(size_t k=0;k<coors.size();k++){
            if(coors.size()!=defAlgs.size() or coors[k]!=defAlgs[k]){
                Algebra alg(defAlgs[k]);
                if(defAlgs[k]=="Q")
                    _iSingle=k;
                else if(alg.isAlgebraOfConsts()){
                    _q.resize(coors.size());
                    _q[k]=alg.val(0.);
                }
                else
                    ABORT(Sstr+"arguments in "+Definition+" must be "+Coors+"or Q and constant, got: "+defAlgs);
            }
        }
    }
    vector<string> var;
    if(Definition.rfind(":")!=string::npos)var=tools::splitString(Definition.substr(Definition.find(":")),':');
    if(var.size()==0)for(size_t k=0;k<defAlgs.size();k++)var.push_back("Q"+tools::str(k));
    if(var.size()!=defAlgs.size())
        ABORT("the number of variables does not match number of arguments in "+Definition);

    for(size_t k=0;k<defAlgs.size();k++){
        defAlgs[k]=tools::substringReplaceAll(defAlgs[k],var[k],"Q");
        _arg.push_back(new Algebra(defAlgs[k]));
        if(_arg.back()==0 or not _arg.back()->isAlgebra())
            ABORT("failed to evaluate algebra "+defAlgs[k]+" as part of "+Definition+"\nFailure: "+Algebra::failures);
    }
}

void AlgebraMulti::checkArgs(const std::vector<std::complex<double> > Q) const {
    if(Q.size()!=_arg.size())ABORT(Str("need")+_arg.size()+"arguments for"+definition()+"got"+Q.size());
}

AlgebraExternal::AlgebraExternal(std::string Definition):AlgebraMulti(Definition,"Eta,Rn"){
    vector<string> par=tools::splitString(tools::stringInBetween(Definition,"[","]"),',');
    if(par.size()!=2)ABORT("need format External[dist,screen](Eta,Rn), got: "+Definition);
    dist=tools::string_to_double(par[0]);
    screen=tools::string_to_double(par[1]);
    trunc.reset(new Algebra("trunc[10,20]"));
}

std::complex<double> AlgebraExternal::valMulti(const std::vector<std::complex<double> > Q) const{
    checkArgs(Q);
    if(abs(Q[1])<3.)return 0;
    complex<double> z=Q[0]*Q[1];
    complex<double> rhoSq=Q[1]*Q[1]-z*z;

    return -trunc->val(Q[1])*(1./sqrt(pow(z-dist,2)+rhoSq+screen));
}

AlgebraCO2Pot::AlgebraCO2Pot(std::string Definition)
    :AlgebraMulti(Definition,"Eta,Rn"),bCO(2.197)
{
    vector<string> par=tools::splitString(tools::stringInBetween(Definition,"[","]"),',');
    if(par.size()!=4)ABORT("need format CO2Pot[runcR,alfa,screenC,screenO](Eta,Rn), got: "+Definition);

    double truncR=Algebra::getParameter(0,Definition);
    trunc.reset(new Algebra("trunc["+tools::str(truncR-5.,2)+","+tools::str(truncR,2)+"]"));
    alfaCharge=Algebra::getParameter(1,Definition);
    screenC=Algebra::getParameter(2,Definition);
    screenO=Algebra::getParameter(3,Definition);
}

std::complex<double> AlgebraCO2Pot::valMulti(const std::vector<std::complex<double> > Q) const{
    checkArgs(Q);
    complex<double>v;
    v =-alfaCharge*(1.+5.*exp(-Q[1]/screenC))/Q[1];
    complex<double> rOff;
    rOff=sqrt(Q[1]*Q[1]+2.*Q[0]*Q[1]*bCO+bCO*bCO);v+=(0.5*(alfaCharge-1.))*(1.+7.*exp(-rOff/screenO))/rOff;
    rOff=sqrt(Q[1]*Q[1]-2.*Q[0]*Q[1]*bCO+bCO*bCO);v+=(0.5*(alfaCharge-1.))*(1.+7.*exp(-rOff/screenO))/rOff;
    v*=trunc->val(Q[1]);
    return v;
}

AlgebraNumerical::AlgebraNumerical(std::string Definition)
    :AlgebraMulti(Definition,"*")
{
    std::string path=tools::stringInBetween(Definition,"Numerical[","]",true);

    std::vector<std::string> mess=
    {std::string(" test potential can be machine generated by creating an input file")
    ,std::string(" at \""+path+"\" with a single line of the format")
    ,std::string(" #   ExternalHO(Eta,Rn):Eta:Rn|[-1,1,21]|[0,10,101]")
    ,std::string(" where ExternalHO(Eta,Rn) is some allowed two-argument algebra")
    ,std::string(" 21 grid points on [-1,1] for Eta")
    ,std::string(" 101 grid points on [0,10] for Rn")
     ,std::string("")
    };

    if(not folder::exists(path)){
        PrintOutput::warning("Potential file \""+path+"\" not found");
        for(auto c: mess)PrintOutput::subTitle(c);
        ABORT("Numerical potential file "+path+" not found");
    }
    AsciiFile file(path);
    std::vector<std::string> comm;
    file.readComments(comm);
    if(comm.size()==1 and comm.front().find("|")!=std::string::npos){
        PrintOutput::message("generate "+Definition+" with "+comm[0]);
        std::vector<std::string> parts=tools::splitString(comm[0],'|');
        std::shared_ptr<const AlgebraMulti> alg(AlgebraMulti::factory(tools::cropString(parts[0])));
        for(size_t k=1;k<parts.size();k++){
            _axes.push_back(tools::rangeToGrid(parts[k]));
            comm.push_back("axis:");
            for(auto g: _axes.back())comm.back()+=" "+(tools::str(g))+",";
            comm.back().pop_back();
        }
        if(_axes.size()!=2)DEVABORT("only for two axes, got: "+comm[0]);
        for(auto g1: _axes[1])
            for(auto g0: _axes[0])
                _values.push_back(alg->valMulti({g0,g1}).real());

        comm.insert(comm.begin()+1,
                    {std::string(" --- machine generated from file with single input line above (see AlgebraMulti.cpp for details ) ---")
                     ,std::string("address as Numerical[full-file-path]]")
                     ,std::string("function values will be linearly interpolated from values here")
                     ,std::string("for points outside axes, 0 is returned")
                     ,std::string("column of values, index of first axis is first (column-wise storage)")
                     ,std::string("format: rows starting with '# axis:' define axes, as below")
                    });
        file.writeComments(comm);
        file.writeCols({_values});
        _axes.clear();
        _values.clear();

    }

    for(auto c: comm){
        if(tools::cropString(c).find("axis:")!=0)continue;
        auto sax=tools::splitString(c.substr(c.find(":")+1),',');
        _axes.push_back(std::vector<double>());
        for(auto v: sax)_axes.back().push_back(tools::string_to_double(v));
    }
    std::vector<std::vector<double>> cols;
    file.readCols(cols);
    if(_axes.size()>2){
        for(auto c: mess)PrintOutput::message(c);
        ABORT("for now, at most 2-dimensional numerical potentials");
    }
    size_t siz(1);
    for(auto a:_axes)siz*=a.size();
    _values=cols[0];
    if(_values.size()!=siz){
        for(auto c: mess)PrintOutput::message(c);
        ABORT(Sstr+"product axis sizes "+siz+"does not match number of values"+_values.size());
    }
}

static double linInt(const std::vector<std::complex<double>> &Q, const std::vector<std::vector<double>> &Ax, const std::vector<double> &Val){

        size_t i=std::upper_bound(Ax.back().begin(),Ax.back().end(),Q.back().real())-Ax.back().begin();

        if(i==0 or i==Ax.back().size())return 0;
        if(not (Ax.back()[i-1]<=Q.back().real() and Q.back().real()<=Ax.back()[i]))ABORT("incorrectly located in axis");

        double v0,v1;
        if(Q.size()==1){
            v0=Val[i-1];
            v1=Val[i];
        }
        else{
            size_t siz(1);
            for(size_t k=0;k<Ax.size()-1;k++)siz*=Ax[k].size();
            v0=linInt({Q.begin(),Q.end()-1},{Ax.begin(),Ax.end()-1},{Val.begin()+((i-1)*siz),Val.end()});
            v1=linInt({Q.begin(),Q.end()-1},{Ax.begin(),Ax.end()-1},{Val.begin()+((i  )*siz),Val.end()});
        }

        return v0+(v1-v0)*(Q.back().real()-Ax.back()[i-1])/(Ax.back()[i]-Ax.back()[i-1]);
}

std::complex<double> AlgebraNumerical::valMulti(const std::vector<std::complex<double> > Q) const{
    if(Q.size()!=_axes.size())ABORT(Sstr+_definition+"needs"+_axes.size()+"arguments, got: "+Q.size());
    std::vector<size_t> idx;
    for(size_t i=0;i<Q.size();i++){
        idx.push_back(std::upper_bound(_axes[i].begin(),_axes[i].end(),Q[i].real())-_axes[i].begin());
        if(idx.back()==0 or idx.back()==_axes[i].size())return 0.;
        if(Q[i].imag()!=0)ABORT(Sstr+"cannot have complex argument"+Q+"at non-zero range of numerical potential"+_definition);
    }
    return linInt(Q,_axes,_values);
}

AlgebraCO2Nuclear::AlgebraCO2Nuclear(std::string Definition)
    :AlgebraMulti(Definition,"Eta,Rn"),bCO(2.197)
{
    if(tools::stringInBetween(Definition,"[","]")!=Definition){
        vector<string> par=tools::splitString(tools::stringInBetween(Definition,"[","]"),',');
        if(par.size()<1)ABORT("need format CO2Nuclear(Eta,Rn) or CO2Nuclear[bond](Eta,Rn) , got: "+Definition);
        bCO=Algebra::getParameter(0,Definition);
    }
}

std::complex<double> AlgebraCO2Nuclear::valMulti(const std::vector<std::complex<double> > Q) const{
    checkArgs(Q);
    complex<double>v(0.);
    if(std::norm(Q[0])>1.)DEVABORT("|Eta|>1");
    complex<double> etaRn2=2.*Q[0]*Q[1];
    v=-6./Q[1]-8./sqrt(Q[1]*Q[1]+bCO*(bCO-etaRn2))-8./sqrt(Q[1]*Q[1]+bCO*(bCO+etaRn2));
    return v;
}

