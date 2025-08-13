// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "potSolid.h"
#include "readInput.h"
#include "constants.h"
#include <map>

using namespace std;

static map<string,const Algebra*> singleSite;

const Algebra* PotSolid::potSingFactory(string Func, int Der){

    string func=Func;
    if(singleSite.count(func+tools::str(Der))==0){
        string para=tools::stringInBetween(Func,"[","]");
        if(para==Func)para="";
        else {
            Algebra c(para);
            if(not c.isAlgebra())ABORT("parameter is not Algebra in "+Func);
            para=tools::str(c.val(0.).real());
        }

        // a few convenient shortcuts
        if(Func=="zero")func="0.";
        if(Func=="sine")func="0.5*(1-cos(pi*Q)";
        if(Func.find("gauss[")==0) func="(1.-exp(-pow[2](Q/"+para+")))";
        if(Func.find("expCos[")==0)func="(1.-exp(-(1-cos(pi*Q))/"+para+"))";

        // norm
        Algebra nrm(func);
        if(not nrm.isAlgebra())ABORT("not an Algebra: "+func);

        string der=func;
        if(Der==1){
            if(Func=="zero")der="0.";
            else if(Func=="sine")der="0.5*pi*sin(pi*Q)";
            else if(Func.find("gauss[")==0) der="(exp(-pow[2](Q/"+para+")))*2*(Q/pow[2]("+para+"))";
            else if(Func.find("expCos[")==0)der="(-exp(-(1-cos(pi*Q))/"+para+"))*(-pi*sin(pi*Q)/"+para+")";
            else ABORT("need to define derivative for "+Func);
        }
        if(nrm.val(1.)!=0.)func="("+der+")*"+tools::str(1./nrm.val(1.).real());
        singleSite[Func+tools::str(Der)]=new Algebra(func);
        if(not singleSite[Func+tools::str(Der)]->isAlgebra())
            ABORT("ill-defined algebra "+func+" failure: "+singleSite[Func+tools::str(Der)]->failures);
    }
    return singleSite[Func+tools::str(Der)];
}

const Algebra* PotSolid::exportFactory(std::string Term){
    std::string file=tools::stringInBetween(Term,"[","]");
    if(Term==file)
        return new PotSolid(ReadInput::main);
    else {
        ReadInput inp(file);
        return new PotSolid(inp);
    }
}
const Algebra* PotSolid::exportDerivativeFactory(std::string Term){
    std::string file=tools::stringInBetween(Term,"[","]");
    if(Term==file)
        return new PotSolid(ReadInput::main,1);
    else {
        ReadInput inp(file);
        return new PotSolid(inp,1);
    }
}


PotSolid::PotSolid(ReadInput & Inp, int Derivative)
    :_der(Derivative)
{
    if(not Inp.found("PotSolid"))return;

    definition()="PotSolid";
    if(Derivative==1)_definition="Derivative_PotSolid";

    if(&Inp!=&ReadInput::main)_definition+="["+Inp.file()+"]";


    //PotSolid: sites, aLattice, height, shape, zeroLevel
    int line=0,nSite=1;
    string shape;
    // the -1 section:
    lim.push_back(0.);     // upper end
    lheight.push_back(0.); // left height
    rheight.push_back(0.); // right height
    zero.push_back(0.);    // zero level
    potSing.push_back(potSingFactory("zero"));
    derSing.push_back(potSingFactory("zero"));

    double heig0=0.,zero0=0.,a;
    string join;
    while(nSite!=0){
        line++;
        Inp.read("PotSolid","sites",nSite,"0","number of sites",line);
        if(line>1 and nSite==0)break;
        Inp.read("PotSolid","aLattice",a,"0.","lattice constant",line);
        Inp.read("PotSolid","shape",shape,"zero",
                 "single site potential: zero,sine,gauss[width],expCos[width]...legacy",line);
        Inp.read("PotSolid","height",heig0,"0.","scale height",line);
        Inp.read("PotSolid","zeroLevel",zero0,"0.","zero level of potential",line);
        Inp.read("PotSolid","join",join,"left","height at joint: left...previous site, right...next site",line);

        for(int k=0;k<nSite;k++){
            if(a==0.)ABORT("cannot have zero lattice constant");
            lim.push_back(lim.back()+a);
            potSing.push_back(potSingFactory(shape));
            derSing.push_back(potSingFactory(shape,1));
            lheight.push_back(heig0); // scale from function value
            rheight.push_back(heig0); // scale from function value
            if(join=="left")      lheight.back()=rheight[rheight.size()-2]-zero.back()+zero0;
            else if(join=="right")rheight[rheight.size()-2]=lheight.back()+zero.back()-zero0;
            else ABORT("illegal value of join="+join);

            zero.push_back(zero0);
        }
    }

    // adjust to zero level at end
    rheight.back()=zero.back();
    double shift=0.5*(lim[0]+lim.back());
    for(int k=0;k<lim.size();k++)lim[k]-=shift;

    if(Derivative==0 and not folder::exists(Inp.outputTopDir()+"potSolid")){
        plot(Inp.outputTopDir()+"potSolid",2*lim[0]-lim[1],lim.back()+lim[1]-lim[0],2001);// oscillatory - take plenty points
        PrintOutput::DEVmessage("on file "+definition());
    }
    if(Derivative==1 and not folder::exists(Inp.outputTopDir()+"derSolid")){
        plot(Inp.outputTopDir()+"derSolid",2*lim[0]-lim[1],lim.back()+lim[1]-lim[0],2001);// oscillatory - take plenty points
        PrintOutput::DEVmessage("on file "+definition());
    }

}

std::complex<double> PotSolid::val(const std::complex<double> Q) const{
    double q=Q.real();
    if(lim[0]<q and q<lim.back() and Q.imag()!=0.)
        ABORT(Str("potential must not be in complex region")+lim[0]+"<="+q+Q+"<="+lim.back());

    if(q<=lim[0])    q=lim[0]    +1.e-11*abs(lim[0]);
    if(q>=lim.back())q=lim.back()-1.e-11*abs(lim.back());

    int site=upper_bound(lim.begin(),lim.end(),q)-lim.begin(); // locate
    double scale=2/(lim[site]-lim[site-1]);
    q=(q-0.5*(lim[site]+lim[site-1]))*scale;     // map into [-1,1]

    if(_der==0){
        if(q<0.)return potSing[site]->val(-q).real()*lheight[site]-zero[site];
        else    return potSing[site]->val( q).real()*rheight[site]-zero[site];
    }
    else if(_der==1){
        if(q<0.)return -derSing[site]->val(-q).real()*lheight[site]*scale;
        else    return  derSing[site]->val( q).real()*rheight[site]*scale;
    }
    else ABORT(Str("at most derivative =1, is ")+_der);
}
