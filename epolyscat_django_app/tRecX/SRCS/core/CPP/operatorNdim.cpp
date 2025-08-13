// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorNdim.h"

#include "stringTools.h"
#include "str.h"

#include "operatorDefinition.h"
#include "operatorData.h"
#include "operatorData.h"
#include "algebra.h"
#include "parameters.h"

using namespace std;

std::map<std::string,OperatorNdim::Kernel> OperatorNdim::special;

OperatorNdim::OperatorNdim(const string &Def, const string &ICoor, const string &JCoor)
{
    setDefaultKernels();
    string def(Def);

    // decompose into single terms
    vector<string> terms=OperatorData::singleTerms(def,ICoor,JCoor);
    complex<double> mult(1.);
    for(size_t k=0;k<terms.size();k++){
        _terms.push_back(Term(terms[k],mult));
    }
}

OperatorNdim::Term::Term(const std::string Single, std::complex<double> Multiplier)
    :_ivd(0),_jvd(0),_mult(Multiplier),_potNdim(0){

    string s(Single);

    size_t lbr=min(s.find("[["),s.find("<"));
    if(lbr==string::npos)ABORT("neither standard nor special operator ");

    if(lbr!=0){
        if(Parameters::isFunction(s.substr(0,lbr)))ABORT("no functions admitted");
        _mult*=*Parameters::pointer(s.substr(0,lbr));
        s=s.substr(lbr);
        lbr=0;
    }

    if(s.find("[[")!=string::npos){
        if(special.count(s)!=0)
            _potNdim=special[s];
        else
            ABORT("special potential not defined: \""+s
                  +"\"\nAvailable potentials:\n"+tools::listMapKeys(special,"\n"));
    }

    else {
        int ncoo=0;
        while (lbr!=string::npos){
            ncoo++;
            string op=s.substr(lbr+1,s.find(">",lbr)-lbr-1);
            lbr=s.find("<",lbr+1);
            if(op.find("d_")==0){
                if(_ivd>0)ABORT("higher derivatives not allowed: "+s);
                _ivd=ncoo;
                op=op.substr(2);
            }
            if(op.length()>2 and op.find("_d")==op.length()-2){
                if(_jvd>0)ABORT("higher derivatives not allowed: "+s);
                _jvd=ncoo;
                op=op.substr(0,op.length()-2);
            }
            if(op=="1" or op=="{}")alg.push_back(0);
            else                   alg.push_back(new Algebra(op));

            // multi-coordinate functions
            if(alg.back()!=0 and not alg.back()->isAlgebra()){
                ABORT("operator not implemented yet: "+Single);
            }
        }
    }
}
std::complex<double>  OperatorNdim::Term::kernel(std::vector<std::complex<double> > Point) const {
    std::complex<double> res=_mult;
    if(_potNdim!=0)
        res*=_potNdim(Point);
    else{
        if(alg.size()!=Point.size())DEVABORT(Sstr+"algebra factors"+alg.size()+" do not match point dimension"+Point.size());
        for(size_t k=0;k<Point.size();k++)
            if(alg[k])res*=alg[k]->val(Point[k]);
    }
    return res;
}

std::string OperatorNdim::Term::str() const {
    Str s("","");
    s+=_mult;
    vector<string>ld(3,"<"),rd(3,">");
    if(_ivd>0)ld[_ivd-1]="<d_";
    if(_jvd>0)rd[_jvd-1]="_d>";
    for(size_t k=0;k<alg.size();k++){
        s+=ld[k];
        if(alg[k]==0)s+"1";
        else s+=alg[k]->definition();
        s+=rd[k];
    }
    return std::move(s);
}

void OperatorNdim::addPotNdim(string Name, Kernel Pot){
    Name="[["+Name+"]]";
    if(special[Name]==0)special[Name]=Pot;
    else if(special[Name]!=Pot){
        ABORT("re-definition of potential");
    }
}













