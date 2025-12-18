// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisSetDef.h"

#include "index.h"
#include "basisAbstract.h"
#include "printOutput.h"
#include "basisExpIm.h"
#include <set>

using namespace std;
using namespace tools;

static std::set<std::string> discreteFunctions({"useIndex","CIion","CIneut","Hybrid","vector"});

string BasisSetDef::size() const {
    string dep=tools::stringInBetween(funcs,"{","}");
    if(dep!=funcs)return tools::str(order)+"-"+dep;
    int siz=order;
    if(not discreteFunctions.count(funcs)){
        bool asymp0=BasisFunction::asympZero(funcs);
        if(first and not (asymp0 and scale<0) and coor.zeroLow)siz--;
        if(last  and not (asymp0 and scale>0) and coor.zeroUp )siz--;
    }
    return tools::str(siz);
}

double BasisSetDef::upBound() const {
    double res=shift+scale;
    if(funcs.find("Channel")==0)return shift+order-1;
    if(discreteFunctions.count(funcs))return order-1;
    if(funcs=="useIndex" or funcs=="CIion" or funcs=="CIneut" or coor.name()=="Orbital" or coor.name()=="Vec")return order-1;
    if(scale<0)return shift;  // negative range element
    if(funcs.find("jacobi")==std::string::npos){
        // hacky mess...
        if(BasisFunction::asympZero(funcs))return DBL_MAX;  // infinit element
    }
    return res;     // standard case
}

double BasisSetDef::lowBound() const {
    if(discreteFunctions.count(funcs))return 0;
    if(scale>0)return shift; // standard case
    // negative range
    if(BasisFunction::asympZero(funcs))return -DBL_MAX; // infinite element
    return shift+scale; // finite negative range element
}


std::string BasisSetDef::str() const{
    // basis definition in a string
    std::string s=name()+" ["+tools::str(lowBound(),3)+","+tools::str(upBound(),3)+"] ("+tools::str(shift,3)+","+tools::str(scale,3)+") "
            +tools::str(order)+"["+size()+"]";
    if(comSca.etaX(shift+scale*0.5)!=1.)s+=tools::str(comSca.etaX(shift+scale*0.5),3);
    s+=" "+modify();
    return s;
}

bool BasisSetDef::operator==(const BasisSetDef &o) const{// comparison: place cheapest and most likely differences first
    if(coor.cString!=o.coor.cString)goto FailComp;
    if(exactIntegral!=o.exactIntegral)goto FailComp;
    if(shift!=o.shift)goto FailComp;
    if(scale!=o.scale)goto FailComp;
    if(par!=o.par)goto FailComp;
    if(first!=o.first)goto FailComp;
    if(last!=o.last)goto FailComp;
    if(not (comSca==o.comSca))goto FailComp;
    if(funcs!=o.funcs)goto FailComp;
    if(margin!=o.margin)goto FailComp;
    if(deriv!=o.deriv)goto FailComp;
    if(order!=o.order)goto FailComp;
    if(minOrder!=o.minOrder)goto FailComp;
    return true;
FailComp:
    if(funcs=="useIndex" and o.funcs=="useIndex")return order==o.order;
    return false;
}

// OBSOLESCENT: do not add new functionality, move to Index::resolveDependence rather than fixing
BasisSetDef BasisSetDef::resolveDependence(const vector<unsigned int> &Pos,const std::vector<const Index *> & Path) const {

    BasisSetDef newDef(*this);
    // curly brackets in square-bracked specification will be ignored
    if(tools::findFirstOutsideBrackets(funcs,"{","[","]")==string::npos)return newDef;

    string dep=tools::stringInBetween(funcs,"{","}");
    // series of special cases
    if(funcs.find("assocLegendre")==0){
        // determine (m) for associated Legendre functions

        if(funcs!=dep){
            //explict constraints:
            // admissible constraints
            vector<string>deps=tools::splitString(dep,'.');
            vector<string> etaDeps=tools::splitString(string("Phi.L-M<.Eta.Lshape=.L-|M|"),'.');
            for(unsigned int k=0;k<deps.size();k++){
                unsigned int l;
                for(l=0;l<etaDeps.size();l++)
                    if(deps[k].find(etaDeps[l])==0)break;
                if(l==etaDeps.size())ABORT("undefined dependency \""+deps[k]
                                           +"\" in funcs="+funcs+"\nadmissible: "+tools::str(etaDeps,3," and "));
            }
        }

        // get the matching axis name
        string phiAx=Path.back()->axisName();
        if(phiAx.find("Eta")==string::npos)ABORT("assocLegendre must be used on Eta-axis, is: "+phiAx);
        phiAx="Phi"+phiAx.substr(3);
        if(dep!=funcs and
                dep.find(phiAx)==string::npos and
                dep.find("Lshape")==string::npos)
            ABORT("axis dependency inconsistent with axis name: "+funcs+" --- no need to specify");

        // locate basis set for axis
        unsigned int k;
        for(k=0;k<Path.size();k++)if(Path[k]->axisName()==phiAx)break;
        if(k==Path.size())ABORT("found not match for "+phiAx+" in axis names");

//        const BasisSet* b=0;
        const BasisTrigon* bExp=dynamic_cast<const BasisTrigon*>(Path[k]->basis());
        if(bExp==0){
            DEVABORT("cannot use BasisSet any longer");
//            b=Path[k]->basisSet();
//            if(b->def.funcs.find("expIm")==string::npos and b->def.funcs.find("cosSin")==string::npos)
//                ABORT("must have expIm or cosSin to match assocLegrendre, found: "+b->str());
        }
        string funcName=bExp->name();

        // determine branch and assign parameters
        unsigned int n=Pos[k];
        int funcPar=bExp->mValueAbs(n);

        // Old sign convention: P_lm = P_l-m
        // is NOT consistent with Gaunt
        // newDef.par.assign(1,abs(b->fs->par[n]));

        // New sign convention: P_lm = (-)^m P_l-m
        //        if(bExp)newDef.par.assign(1,std::abs(bExp->mValue(n)));
        //        else    newDef.par.assign(1,b->fs->par[n]);
        newDef.par.assign(1,std::abs(funcPar));
        newDef.order-=min(order,(unsigned int)abs(newDef.par[0]));

        // here we may add a lower constraint as for L-M<....
        if(dep.find("L-M<")!=string::npos){
            PrintOutput::message("Deprecated old constraints. Use IndexConstraint instead");
            //            if(b->def.funcs.find("expIm")==string::npos)
            if(funcName.find("expIm")==string::npos)
                ABORT("L-M constraint only for expIm functions on Phi-coordinate, is: "+funcName);
            string lm=dep.substr(dep.find("L-M<")+4);
            int diff;
            diff=tools::string_to_int(lm.substr(0,lm.find(".")));
            newDef.order=min(newDef.order,(unsigned int)max(1,diff-abs(int(newDef.par[0]))+funcPar));
        }
        if(dep.find("L-|M|<")!=string::npos){
            PrintOutput::message("Deprecated old constraints. Use IndexConstraint instead");
            if(funcName.find("expIm")==string::npos)
                ABORT("L-|M| constraint only for expIm functions on Phi-coordinate, is: "+funcName);
            string lm=dep.substr(dep.find("L-|M|<")+6);
            int diff=tools::string_to_int(lm.substr(0,lm.find(".")));
            newDef.order=min(newDef.order,(unsigned int)max(1,diff));
        }
    }

    // no dependency in definition
    if(dep==funcs)return newDef;

    if(funcs.find("expIm")==0)
    {
        // constraint m1+m2=M
        if(dep.find(".M=")==string::npos)ABORT("undefined dependence in: "+funcs+", example format expIm{Phi1.M=0}");
        ABORT("Obsolete. Use BasisConstraint: kind=M[a,b], axes=... instead");

//        // get the matching axis name and basis set
//        unsigned int k;
//        string phiAx=dep.substr(0,dep.find(".M="));
//        for(k=0;k<Path.size();k++)if(Path[k]->axisName()==phiAx)break;
//        if(k+1>Path.size())ABORT("found no match for "+phiAx+" in axis names");
//        const BasisSet* b=Path[k]->basisSet();
//        if(Path[k]->strNode().find("expIm")==string::npos
//                and Path.back()->strNode().find("expIm")==string::npos)
//            ABORT("must have both expIm functions for M=0 constraint");

//        // compute parameter
//        newDef.par.assign(1,tools::string_to_double(dep.substr(dep.find("=")+1))-b->fs->par[Pos[k]]);
//        if(order/2>=(unsigned int)abs(newDef.par[0]))newDef.order=1;
//        else                                         newDef.order=0; // required m outside of range
    }

    if(funcs.find("assocLegendre")==0){
        if(dep.find("Lshape=")!=string::npos)
        {
            DEVABORT("Deprecated old constraints. Use IndexConstraint instead");
//            unsigned int k;
//            for(k=0;k<Path.size();k++)if(Path[k]->axisName()=="Eta1")break;
//            if(k+1>Path.size())ABORT("found not Eta1 above curren, required for Lshape=");
//            const BasisSet* b=Path[k]->basisSet();
//            if(b->def.funcs.find("assocLegendre")==string::npos and Path.back()->basisSet()->def.funcs.find("assocLegendre")==string::npos)
//                ABORT("must have both expIm functions for 'Lshape=' constraint");
//            unsigned int l1=abs(b->param(0))+Pos[k];
//            int lWidth=tools::string_to_int(dep.substr(dep.find("=")+1));
//            if(lWidth<1)ABORT("must admit at least one angular function, Lshape>0, is: "+funcs);
//            if(l1>=lWidth)newDef.order=max(lWidth-1,int(newDef.par[0])+1);
        }
        else
            if(dep.find("Phi")==string::npos)
                ABORT("undefined dependence in: "+funcs+", example format assocLegendre{Lshape=5}");


    }
    else if(funcs.find("besselCoulomb")==0){
        DEVABORT("besselCoulomb no longer available");
//        // get the matching axis name
//        string etaAx=Path.back()->axisName();
//        if(etaAx.find("Rn")==string::npos)ABORT("besselCoulomb must be used on Rn-axis, is: "+etaAx);
//        etaAx="Eta"+etaAx.substr(2);

//        // locate basis set for axis
//        unsigned int k;
//        for(k=0;k<Path.size();k++)if(Path[k]->axisName()==etaAx)break;
//        if(k==Path.size())ABORT("found no match for "+etaAx+" in axis names");
//        const BasisSet* b=Path[k]->basisSet();

//        // determine branch and assign parameters
//        unsigned int n=Pos[k];

//        // assign angular momentum to par
//        newDef.par.assign(1,b->physical(n));
    }

    return newDef;
}
