// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "operatorDefinition.h"

#include "tools.h"
#include "readInput.h"
#include "algebra.h"
#include "index.h"
#include "str.h"
#include "parameters.h"
#include "basisAbstract.h"

#include "operatorFloorEE.h"
using namespace std;

bool OperatorDefinition::skipAxis(const Index* Idx){
    // these axis do not carry operators
    if(Idx->depthOfDuplicate()!=Index::npos)return true;
    if(Idx->axisName().find("surf")==0)return true;
    if(Idx->continuity()!=Index::npos && Idx->axisName().find("Ion")!=0)return true;
    if(Idx->axisName().find("spec")==0){
        for(const Index* idx=Idx;idx!=0;idx=idx->descend())
            if("k"+Idx->axisName().substr(4)==idx->axisName())return true;
    }
    return false;
}

const OperatorDefinition OperatorDefinition::identity = OperatorDefinition("[[Identity]]");

std::map<std::string,std::map<std::string,std::string> > OperatorDefinition::standardOperators;

OperatorDefinition &OperatorDefinition::construct(string Def, string Hierarchy){
    // replaces old constructor...

    OperatorDefinition inDef(Def);
    if(Def.find("<<1>>")!=string::npos)ABORT("OBSOLETE operator defintion <<1>>, replace by <<Overlap>>");

    string coors="",def="undefined";
    //    vector<string> defs=OperatorData::terms(Def);
    std::vector<OperatorDefinition> defs=inDef.terms();
    if(tools::subStringCount(Def,"<<")>1){
        string s;
        for(size_t k=0;k<defs.size();k++){
            if(defs[k].find("<<")==string::npos)s+=defs[k];
            else s+=(OperatorDefinition().construct(defs[k],Hierarchy));
        }
        assign(s);
        return *this;
    }

    // try to resolve pre-defined operator
    if(Hierarchy!=""){

        coors=Index::coordinates(Hierarchy);
        def=OperatorDefinition::get(Def,coors);
        // remove constraints
        if(Def.find("$")!=std::string::npos)
            PrintOutput::DEVwarning("constraints removed - may dramatically increase setup time");
        Def=Def.substr(0,def.find("$"));
    }

    // take literally
    if(def.find("undefined")==string::npos){
        assign(def);
        return *this;
    }

    // if Hierachy contains FE-axis, move factors properly
    if(Hierarchy!="" and coors!=Hierarchy){
        //        vector<string> term=OperatorData::terms(Def);
        std::vector<OperatorDefinition> term=inDef.terms();
        string s;
        for(size_t k=0;k<term.size();k++){
            s+=OperatorDefinition(term[k]).axisToCoor(Hierarchy);
        }
        assign(s);
    }
    else
        assign(Def);
    return *this;
}

OperatorDefinition & OperatorDefinition::tsurffDropTerms(const Index *IIndex, const Index *JIndex){
    std::vector<OperatorDefinition> trms=terms();
    if(IIndex->hierarchy().find("spec")!=string::npos or JIndex->hierarchy().find("spec")!=string::npos){
        dropTerms(trms,IIndex->hierarchy());
    }
    clear();
    for(unsigned int k=0;k<trms.size();k++)*this+=trms[k];
    return *this;
}

void OperatorDefinition::dropTerms(std::vector<string> &Terms, const Index* IIndex, const Index *JIndex)
{
    if(IIndex==0)return;
    if(IIndex!=JIndex)DEVABORT("dropTerms only for operators mapping space onto itself");
    vector<string> terms1(0);
    std::vector<OperatorDefinition> terms;
    for(auto t: Terms)terms.push_back(OperatorDefinition(t));

    for(unsigned int k=0;k<terms.size();k++){
        // Go through all terms
        bool drop = false;
        string temp = terms[k];

        const Index *iI=IIndex,*jI=JIndex;
        while(not drop and not iI->isLeaf() and not jI->isLeaf()){
            if(iI->depthOfDuplicate()!=Index::npos){
                // Continuity level - do nothing as there is no string associated.
                if(iI->depthOfDuplicate()!=jI->depthOfDuplicate())ABORT("continuity levels do not match");
            }
            else if((iI->axisName().substr(0,4)=="surf" and jI->axisName().substr(0,4)=="surf") or
                    (iI->axisName().substr(0,4)=="spec" and jI->axisName().substr(0,4)=="spec")){
                // Surface/Spec level - do nothing as there is no string associated.
            }
            else {
                //                string a = OperatorData::first(terms[k],true);
                string a = terms[k].first(true);
                terms[k]=remainder(terms[k]);
                if(iI->axisName().substr(0,1)=="k" and jI->axisName().substr(0,1)=="k"){ // Momentum index
                    if(a.find("<1>")==string::npos) {drop=true; break;}
                }
            }
            iI=iI->child(0);
            jI=jI->child(0);
        }
        if(not drop)terms1.push_back(temp);
    }

    terms.clear();
    for(auto t: terms1)terms.push_back(OperatorDefinition(t));
    //    terms = terms1;
}
void OperatorDefinition::dropTerms(std::vector<string> &terms, std::string Hierarchy){
    std::vector<OperatorDefinition> trms;
    for(auto t: terms)trms.push_back(OperatorDefinition(t));
    dropTerms(trms,Hierarchy);
    terms.clear();
    for(auto t: trms)terms.push_back(t);
}


static bool multiElectron(std::string Hierarchy){
    if(Hierarchy.find("Ion")!=std::string::npos)return true;
    if(Hierarchy.find("Neutral")!=std::string::npos)return true;
    if(tools::subStringCount(Index::coordinates(Hierarchy),"Rn")>1)return true;
    if(tools::subStringCount(Index::coordinates(Hierarchy),"Eta")>1)return true;
    return false;
}

void OperatorDefinition::dropTerms(std::vector<OperatorDefinition> &terms, std::string Hierarchy)
{
    vector<OperatorDefinition> terms1;
    for(OperatorDefinition termK: terms){
        // decide acceptance of term
        bool accept;
        if((termK.find("EEInteraction")!=std::string::npos or termK.find("eeInt6DHelium")!=std::string::npos)
                and not multiElectron(Hierarchy))
            accept=false;
        else if(not termK.isExpanded() or not termK.isStandard(0,0)){
            // special operators handle k- and spec-factors only if there is matching <1> or <Id>

            // leading standard factors (if any)
            std::vector<OperatorDefinition>fac;
            OperatorDefinition remain(termK);
            for(fac.push_back(remain.first());fac.back()[0]=='<' and fac.back()[1]!='<';remain=remainder(remain));

            std::vector<std::string> coors=tools::splitString(Index::coordinates(Hierarchy),'.');
            accept=true;
            for(size_t k=0;k<coors.size() and accept ;k++){
                if(coors[k][0]=='k' or coors[k].find("spec")==0){
                    accept=fac.size()>k and (fac[k]=="<1>" or fac[k]=="<Id>" or fac[k].find("<GridWeight")!=std::string::npos);
                }
            }
        }
        else {
            // keep only terms where k- and spec-factors are <1> or <Id>
            std::string coor=Index::coordinates(Hierarchy)+"."; // reduce hierarchy to actual coordinates
            // check whether factor count matches coordinates
            if(std::count(coor.begin(),coor.end(),'.')!=std::count(termK.begin(),termK.end(),'<')
                    and termK.find("[[")==string::npos and termK.find("<<")==string::npos){
                DEVABORT(termK+" does not match: "+coor+": need as many <..>'s as there are coordinates");
            }

            // k-factor and spec-factors need <1> or <Id>
            for(size_t pos0=termK.find("<");coor!="" and pos0<termK.find("<<");coor=coor.substr(coor.find(".")+1),pos0=termK.find("<",pos0+2)){
                if(coor[0]=='k' and pos0==string::npos)break; // not tensor product structure - cannot select
                if((coor[0]=='k' or coor.substr(0,4)=="spec") and
                        not(termK.find("<1>",pos0)==pos0
                            or termK.find("<Id>" ,pos0)==pos0
                            or termK.find("<GridWeight",pos0)==pos0)
                        )break;
            }
            accept=coor=="";
        }

        if(accept)terms1.push_back(termK);
    }
    terms=terms1;
}

std::string OperatorDefinition::parentHierarchy(std::string Hierarchy){
    // given a hierarchy, convert special names to where they derived from
    // note: historically, there is some inconsistency in the use of surf/ValDer and spec/k
    //       this is fixed here.
    Hierarchy+=".";
    string parent;
    for(size_t i=0,j=Hierarchy.find('.');j!=string::npos;i=j+1,j=Hierarchy.find('.',i)){
        if     (Hierarchy.find("spec",i)==i)parent+=Hierarchy.substr(i+4,j-i-3);
        else if(Hierarchy.find("surf",i)==i){
            string ax=Hierarchy.substr(i+4,j-i-3);
            if(Hierarchy.find("ValDer"+ax)!=string::npos)parent+="ValDer"+ax;
            else                                         parent+=ax;
        }
        else if(Hierarchy.find("k",i)==i)parent+=Hierarchy.substr(i+1,j-i);
        else parent+=Hierarchy.substr(i,j-i+1);
    }
    return parent.substr(0,parent.length()-1);
}

// this is a rather shaky un-bracketing function
OperatorDefinition OperatorDefinition::unBracket(string Term,string Front,string Back){
    if(Term=="")return OperatorDefinition("");

    // decompose into terms
    vector<std::string> terms,signs,fac,sFac;
    tools::splitString(Term,"+-",terms,signs,"<(",">)");

    tools::splitString(Front,"+-",fac,sFac,"<(",">)");
    if(signs[0]==" ")signs[0]="+";
    if(sFac.size()==0){
        fac.push_back("");
        sFac.push_back("");
    }
    if(sFac[0]=="")sFac[0]="+";

    // remove outermost bracket (if any) from each term
    string res;
    for(unsigned int k=0;k<terms.size();k++){
        // distribute qualifiers
        string back(Back);
        size_t qPos=tools::findLastOutsideBrackets(terms[k],"$","(",")");
        if(qPos!=string::npos){
            if(Back!="")ABORT("multiple qualifiers in Term:\n"+Term);
            back=terms[k].substr(qPos);
            terms[k]=terms[k].substr(0,qPos);
        }
        if(signs[k]==sFac[0])signs[k]="+";
        else                 signs[k]="-";

        size_t i0=tools::findFirstOutsideBrackets(terms[k],"(","<",">");
        if(i0==string::npos)
            res+=signs[k]+fac[0]+terms[k]+back;
        else{
            size_t i1=tools::findLastOutsideBrackets(terms[k],")","<",">");
            if(i1+1<qPos)back=terms[k].substr(i1+1,qPos-i1-1); // new
            string fact=signs[k]+tools::cropString(fac[0])+tools::cropString(terms[k].substr(0,i0));
            string termk=terms[k].substr(i0+1,i1-i0-1);
            res+=unBracket(termk,fact,back);
        }
    }
    return res;
}

/// return definition with factors at FE levels moved to end
OperatorDefinition OperatorDefinition::axisToCoor(const std::string Hierarchy){
    //Warning: remove the setup constraints for now
    *this=substr(0,find("$"));

    if(not isStandard(0,0))return *this; // do not modify special operators
    size_t i0=find("<");
    size_t inext=find("<",i0+1);
    if(inext==string::npos)return *this; // single factor, no rearangement
    if(inext==i0+1)return *this; // predefined operator, no rearrangement

    size_t i1=find(">");
    size_t first=Hierarchy.find(".");
    if(first==string::npos)DEVABORT("hierarchy shorter than operator factors?:\n"+*this);

    if(Hierarchy.find(Hierarchy.substr(0,first),first)==string::npos){
        // not FE level, continue to next factor
        return substr(0,i1+1)+OperatorDefinition(substr(i1+1)).axisToCoor(Hierarchy.substr(first+1)); // not fem axis
    }
    else {
        // FE level: move current factor to end and continue to next
        OperatorDefinition s(*this);
        int lenDef=min(s.length(),s.find("$"));
        s.insert(s.begin()+lenDef,s.begin()+i0,s.begin()+i1+1); // copy to end
        if(i1!=find("><") and i1!=find(">iLaserA")){
            // move factor to front
            string par0=substr(0,find("<"));
            string par1=substr(i1+1,find("<",i1+1)-i1-1);
            ABORT("multiple factors in "+*this+" "+par0+" "+par1);
        }
        s.replace(s.begin()+i0,s.begin()+i1+1,""); // remove from original location
        string constraint;
        return s.substr(0,i0)+OperatorDefinition(s.substr(i0)).axisToCoor(Hierarchy.substr(first+1))+constraint;
    }
}

bool OperatorDefinition::isLocal(const Index *IIndex, const Index *JIndex) const {
    if(not isStandard(IIndex, JIndex)){
        if(*this == "[[eeInt6DHelium]]") return true;
        if(*this == "[[Pot3d]]") return true;
        DEVABORT("define isLocal for non-standard operator "+*this);
    }
    // exchange term is non-local but has standard format
    if(find("XC")!=string::npos)return false;
    return true;
}


bool OperatorDefinition::isStandard(const Index *IIndex, const Index *JIndex) const {
    if(IIndex and IIndex->basis()->ndim())return false;
    if(JIndex and JIndex->basis()->ndim())return false;
    if(find(">")==string::npos) return false;
    vector<string>terms,signs;
    tools::splitString(*this,"+-",terms,signs,"<|",">>");
    return terms.size()!=0 and terms[0].find_first_of("<|")!=string::npos;
}

void OperatorDefinition::specialConstrain(UseMatrix& Mult, const Index *IIndex, const Index *JIndex) const{
    if(*this == "[[eeInt6DHelium]]"){
        OperatorFloorEE::constrain(Mult, IIndex, JIndex);
    }
}

bool OperatorDefinition::isSeparable() const {
    return (tools::subStringCount(*this,":","<([",">)]")!=0);
}
OperatorDefinition OperatorDefinition::remainder(const OperatorDefinition & Def)
{
    if(not Def.isStandard(0,0))return Def;
    return Def.find(">")<Def.find("<<")?OperatorDefinition(Def.substr(Def.find(">")+1)):Def;
}

OperatorDefinition OperatorDefinition::constrain(UseMatrix &Mult, const Index *IIndex, const Index *JIndex) const {
    if(IIndex->basis()->hybrid() and JIndex->basis()->hybrid())return remainder(*this);
    if(IIndex->axisName()!=JIndex->axisName())return *this;
    OperatorDefinition remains(*this);
    size_t itype=find("$");
    int band=INT_MAX;
    if(itype!=string::npos){

        // sanity checks for constraints
        size_t idot0=find(".",itype);
        if(substr(itype,idot0-itype)!="$BAND")ABORT("Constraint not implemented, need \"$BAND\": "+*this);

        // band-width of constraint
        size_t idot1=min(length(),find(".",idot0+1));
        if(idot0<idot1){
            if(substr(idot0+1,idot1-idot0-1)=="*")
                band=INT_MAX;
            else
                band=tools::string_to_int(substr(idot0+1,idot1-idot0-1));
        }

        // remove present constraint from definintion
        remains.erase(remains.begin()+idot0,remains.begin()+idot1);

        // sanity checks for constraints
        string inam=IIndex->basis()->name(),jnam=JIndex->basis()->name();
        if(IIndex->axisName()=="Phi"){
            if((inam.find("PlaneWave")==string::npos and inam.find("CosSin")==string::npos and inam.find("ExpIm")==string::npos) or
                    (jnam.find("PlaneWave")==string::npos and jnam.find("CosSin")==string::npos and jnam.find("ExpIm")==string::npos))
            {
                PrintOutput::warning("Phi-constraint not applied, not trigonometric basis: "+inam+" | "+jnam,5);
                band=INT_MAX;
            }
        }
        if(IIndex->axisName()=="Eta"){
            if(inam.find("assocLegendre")==string::npos or jnam.find("assocLegendre")==string::npos){
                PrintOutput::warning("Eta-constraint not applied, not associatedLegendre basis: "+inam+" | "+jnam,5);
                band=INT_MAX;
            }
        }
    }

    if(band<INT_MAX)
        for(size_t c=0;c<Mult.cols();c++)
            for(size_t r=0;r<Mult.rows();r++)
                if(abs(int(r)-int(c))>band)Mult(r,c)=0.;

    if(not isStandard(IIndex, JIndex)){
        specialConstrain(Mult, IIndex, JIndex);
    }

    return remainder(remains);
}

static string operatorBrackets("([<");
string OperatorDefinition::parameter(OperatorDefinition &Remainder) const{
    size_t bop=find_first_of(operatorBrackets),bra=0;

    while(bra<bop){
        if(bop==string::npos)ABORT("not an operator: "+*this);
        bra=bop;
        switch (string(*this)[bra]){
        case '<': break; // anything that starts with '<' is an operator
        case '[':
            if(find("[[")!=bop)bop=find_first_of(operatorBrackets,bra+1); // only [[ begins an opererator
            break;
        case '(':
            size_t clos=find(')');
            if(clos==string::npos)ABORT("unmatched ')' in "+*this);
            if(min(find('<',bra+1),find("[[",bra+1))>clos)
                bop=find_first_of(operatorBrackets,bra+1);
            break;
        }

        // make sure that we can alias &Remainder with this
        string fac=substr(0,bop);
        Remainder=substr(bop);
        if(fac.find('?')!=std::string::npos){
            fac=tools::splitString(fac,'?')[1];
        }
        return tools::cropString(fac);
    }
    DEVABORT("parameter undefined");
}

void OperatorDefinition::checkVariable(std::string Message,std::string Name) const {
    vector<string> nonConst,timeDep;
    for(unsigned int k=0;k<Parameters::table.size();k++)
        if(Parameters::isFunction(Parameters::table[k].name) and find(Parameters::table[k].name)!=string::npos){
            if(Parameters::table[k].name.find("[t]")!=string::npos)timeDep.push_back(Parameters::table[k].name);
            else                                                  nonConst.push_back(Parameters::table[k].name);
        }
    if(nonConst.size()+timeDep.size()>0){
        PrintOutput::newLine();
        PrintOutput::message("non-constant parameters in operator \""+Name+"\": "+tools::str(timeDep,3,", ")+tools::str(nonConst,3,", "));

        if(nonConst.size()>0 and Message!="")PrintOutput::warning(Message);
    }
}

void OperatorDefinition::setup(){
    if(standardOperators.size()>0)return; // already set up
    standardOperators["GridWeights"];

    standardOperators["Overlap"]["Phi.PXi.PEta"]="0.25<1><Q><1>+0.25<1><1><Q>";
    standardOperators["Laplacian"]["Phi.PXi.PEta"]="<1>(<d_Q_d><1>+<1><d_Q_d>)+0.25<d_1_d>(<1/Q><1>+<1><1/Q>)";
    standardOperators["Coulomb"]["Phi.PXi.PEta"]="0.5<1><1><1>";
    standardOperators["Z"]["Phi.PXi.PEta"]="(<0.125><Q*Q><1>-<0.125><1><Q*Q>)";
    standardOperators["D/DZ"]["Phi.PXi.PEta"]="<0.5>(<Q/2_d><1>-<1><Q/2_d>-<d_Q/2><1>+<1><d_Q/2>)";
    //  this below is wrong
    //    standardOperators["D/DY"]["Phi.PXi.PEta"]=
    //              "(-1/8)(<sqrt(Q)><1/sqrt(Q)>+<1/sqrt(Q)><sqrt(Q)>)(<cos(Q)_d>-<d_cos(Q)>)"
    //              "(1/4)(<sqrt(Q)_d><sqrt(Q)>+<sqrt(Q)><sqrt(Q)_d>-<d_sqrt(Q)><sqrt(Q)>-<sqrt(Q)><d_sqrt(Q)>)<sin(Q)>";

    //    Coriolis term for rotating frame of reference
    standardOperators["CoriolisTheta"]["Phi.PXi.PEta"]=
            "-<0.125*sin(Q)_d><sqrt(Q*Q*Q)><1/sqrt(Q)>+<0.125*sin(Q)_d><1/sqrt(Q)><sqrt(Q*Q*Q)>"
            "-<0.25*cos(Q)><sqrt(Q*Q*Q)_d><sqrt(Q)>-<0.25*cos(Q)><sqrt(Q)_d><sqrt(Q*Q*Q)>"
            "+<0.25*cos(Q)><sqrt(Q*Q*Q)><sqrt(Q)_d>+<0.25*cos(Q)><sqrt(Q)><sqrt(Q*Q*Q)_d>"
            ;

    standardOperators["Id"]["Rn.Phi.Eta"]="<1><1><1>$BAND.0.0.0";
    standardOperators["Id"]["Vec.Ion"]="<1><1>";
    standardOperators["Id"]["Ion"]="<1>";
    standardOperators["Id"]["Rn1.Phi1.Eta1.Rn2.Phi2.Eta2"]="<1><1><1><1><1><1>+<1><1><1><1><1><1>";

    standardOperators["FullOverlap"]["Ion.Rn.Phi.Eta"]="<1><1><1><1>+<allOnes><NONORTH1><allOnes><allOnes>";
    standardOperators["FullOverlap"]["Neutral"]="0.5(<1>+<1>)";
    standardOperators["EEInteraction"]["Ion.Rn.Phi.Eta"]="+<EEInt><1><1><1>$BAND.0.0.0.0+<allOnes><NONORTHEEInt><allOnes><allOnes>";
    standardOperators["TestOp"]["Ion.Rn.Phi.Eta"]="<allOnes><NONORTHLaplacian><allOnes><allOnes>";

    standardOperators["Laplacian"]["Ion.Rn.Phi.Eta"]="<1><d_1_d><1><1>$BAND.0.0.0.0+<1><1/(Q*Q)><1><d_(1-Q*Q)_d>$BAND.0.0.0.0+<1><1/(Q*Q)><d_1_d><1/(1-Q*Q)>$BAND.0.0.0.0"
            +string("+<Laplacian><1><1><1>$BAND.0.0.0.0+<allOnes><NONORTHLaplacian><allOnes><allOnes>");
    standardOperators["Laplacian"]["Ion"]=standardOperators["Laplacian"]["Neutral"]="<Laplacian>";
    standardOperators["Laplacian"]["Vec.Ion"]="<1><Laplacian>";
    standardOperators["Laplacian"]["X"]=standardOperators["Laplacian"]["Y"]=standardOperators["Laplacian"]["Z"]="<d_1_d>";
    standardOperators["Laplacian"]["X.Y"]="<1><d_1_d>+<d_1_d><1>";
    standardOperators["Laplacian"]["X.Y.Z"]="<1><d_1_d><1>+<d_1_d><1><1>+<1><1><d_1_d>";

    standardOperators["EEInteraction"]["Ion"]=standardOperators["EEInteraction"]["Neutral"]="<EEInt>";
    standardOperators["EEInteraction"]["Vec.Ion"]="<1><EEInt>";

    standardOperators["AngularMomentumSquared"]["Rn.Phi.Eta"]="<1><1><d_(1-Q*Q)_d>$BAND.0.0.0+<1><d_1_d><1/(1-Q*Q)>$BAND.0.0.0";
    standardOperators["AngularMomentumZsquared"]["Rn.Phi.Eta"]="<1><d_1_d><1>$BAND.0.0.0";
    standardOperators["Parity"]["Rn.Phi.Eta"]="<1><Parity><Parity>$BAND.0.0.0";

    standardOperators["AngularMomentumZ"]["Phi.Eta"]="<-i_d><1>";
    standardOperators["AngularMomentumX"]["Phi.Eta"]="(<i*cos(Q)_d><Q/sqrt(1-Q*Q)>-<i*sin(Q)><sqrt(1-Q*Q)_d>)";
    standardOperators["AngularMomentumY"]["Phi.Eta"]="(<i*sin(Q)_d><Q/sqrt(1-Q*Q)>+<i*cos(Q)><sqrt(1-Q*Q)_d>)";
    // Coriolis term for rotating frame of reference
    standardOperators["CoriolisTheta"]["Rn.Phi.Eta"]="<i>"+standardOperators["AngularMomentumY"]["Phi.Eta"]+"$BAND.0.4.4";

    // these should be made automatic
    standardOperators["AngularMomentumZ"]["Phi.Eta.Vec"]=standardOperators["AngularMomentumZ"]["Phi.Eta"]+"<Id>$BAND.4.4.0";
    standardOperators["AngularMomentumX"]["Phi.Eta.Vec"]=standardOperators["AngularMomentumX"]["Phi.Eta"]+"<Id>$BAND.4.4.0";
    standardOperators["AngularMomentumY"]["Phi.Eta.Vec"]=standardOperators["AngularMomentumY"]["Phi.Eta"]+"<Id>$BAND.4.4.0";
    // these should be made automatic
    standardOperators["AngularMomentumZ"]["Phi.Eta.Rn"]=standardOperators["AngularMomentumZ"]["Phi.Eta.Vec"];
    standardOperators["AngularMomentumX"]["Phi.Eta.Rn"]=standardOperators["AngularMomentumX"]["Phi.Eta.Vec"];
    standardOperators["AngularMomentumY"]["Phi.Eta.Rn"]=standardOperators["AngularMomentumY"]["Phi.Eta.Vec"];
    standardOperators["Laplacian"]["Phi.Rho"]="<1><d_1_d>$BAND.0.*+<d_1_d><1/(Q*Q)>";

    standardOperators["Laplacian"]["Rn.Phi.Eta"]="<d_1_d><1><1>$BAND.0.0.0+<1/(Q*Q)><1><d_(1-Q*Q)_d>$BAND.0.0.0+<1/(Q*Q)><d_1_d><1/(1-Q*Q)>$BAND.0.0.0";
    standardOperators["Laplacian"]["Rn1.Phi1.Eta1"]=standardOperators["Laplacian"]["Rn2.Phi2.Eta2"]=standardOperators["Laplacian"]["Rn.Phi.Eta"];
    standardOperators["Laplacian"]["Vec.Rn2.Phi2.Eta2"]="<1><d_1_d><1><1>$BAND.0.0.0.0+<1><1/(Q*Q)><1><d_(1-Q*Q)_d>$BAND.0.0.0.0+<1><1/(Q*Q)><d_1_d><1/(1-Q*Q)>$BAND.0.0.0.0";

    standardOperators["Laplacian"]["Rn1.Phi1.Eta1.Rn2.Phi2.Eta2"]=  "<d_1_d><1><1><1><1><1>+<1/(Q*Q)><1><d_(1-Q*Q)_d><1><1><1>+<1/(Q*Q)><d_1_d><1/(1-Q*Q)><1><1><1>";
    standardOperators["Laplacian"]["Rn1.Phi1.Eta1.Rn2.Phi2.Eta2"]+="+<1><1><1><d_1_d><1><1>+<1><1><1><1/(Q*Q)><1><d_(1-Q*Q)_d>+<1><1><1><1/(Q*Q)><d_1_d><1/(1-Q*Q)>";
    standardOperators["Laplacian"]["Rn1.Phi1.Eta1.nSurface.Phi2.Eta2"]=  "<d_1_d><1><1><1><1><1>+<1/(Q*Q)><1><d_(1-Q*Q)_d><1><1><1>+<1/(Q*Q)><d_1_d><1/(1-Q*Q)><1><1><1>";
    standardOperators["Laplacian"]["nSurface.Phi1.Eta1.Rn2.Phi2.Eta2"]=  "<1><1><1><d_1_d><1><1>+<1><1><1><1/(Q*Q)><1><d_(1-Q*Q)_d>+<1><1><1><1/(Q*Q)><d_1_d><1/(1-Q*Q)>";

    standardOperators["Parabolic"]["Rn.Phi.Eta"]="<Q*Q><1><1>";
    standardOperators["Parabolic"]["Rn1.Phi1.Eta1.Rn2.Phi2.Eta2"]="<Q*Q><1><1><1><1><1>+<1><1><1><Q*Q><1><1>";
    standardOperators["Parabolic"]["X"]=standardOperators["Parabolic"]["Y"]=standardOperators["Parabolic"]["Z"]="<Q*Q>";
    standardOperators["Parabolic"]["X.Y"]="<Q*Q><1>+<1><Q*Q>";
    standardOperators["Parabolic"]["X.Y.Z"]="<Q*Q><1><1>+<1><Q*Q><1>+<1><1><Q*Q>";

    string smoothFunction("1"),smoothDerivative("0");
    if(Algebra::isSpecialConstant("Rtrunc")){
        if(not Algebra::isSpecialConstant("Rsmooth"))ABORT("for using potential cutoff radius Rcutoff, must also define Rsmooth");
        smoothFunction  ="trunc[Rsmooth,Rtrunc]";
        smoothDerivative="trunc[Rsmooth,Rtrunc,1]";
    }

    standardOperators["Coulomb"]["Rn1.Phi1.Eta1"]= standardOperators["Coulomb"]["Rn2.Phi2.Eta2"] =standardOperators["Coulomb"]["Rn.Phi.Eta"]="<"+smoothFunction+"/Q><1><1>$BAND.0.0.0";
    standardOperators["Coulomb"]["Vec.Rn1.Phi1.Eta1"]=standardOperators["Coulomb"]["Vec.Rn2.Phi2.Eta2"]="<1><"+smoothFunction+"/Q><1><1>$BAND.0.0.0.0";
    standardOperators["Coulomb"]["Ion.Rn.Phi.Eta"]="-2.<1><"+smoothFunction+"/Q><1><1>$BAND.0.0.0.0+<NucCoulomb><1><1><1>$BAND.0.0.0.0+<allOnes><NONORTHCoulomb><allOnes><allOnes>";
    standardOperators["Coulomb"]["Ion"]=standardOperators["Coulomb"]["Neutral"]="<NucCoulomb>";
    standardOperators["Coulomb"]["Vec.Ion"]="<1><NucCoulomb>";


    standardOperators["Coulomb"]["Rn1.Phi1.Eta1.Rn2.Phi2.Eta2"]="<"+smoothFunction+"/Q><1><1><1><1><1>+<1><1><1><"+smoothFunction+"/Q><1><1>";
    standardOperators["Coulomb"]["Rn1.Phi1.Eta1.nSurface.Phi2.Eta2"]="<"+smoothFunction+"/Q><1><1><1><1><1>";
    standardOperators["Coulomb"]["nSurface.Phi1.Eta1.Rn2.Phi2.Eta2"]="<1><1><1><"+smoothFunction+"/Q><1><1>";

    standardOperators["X"]["Rn.Phi.Eta"]="<Q><cos(Q)><sqrt(1-Q*Q)>";
    standardOperators["Y"]["Rn.Phi.Eta"]="<Q><sin(Q)><sqrt(1-Q*Q)>";
    standardOperators["Z"]["Rn.Phi.Eta"]="<Q><1><Q>";

    // derivative operators in explicitly anti-symmetric form
    string ddxPolar="(<0.5_d><cos(Q)><sqrt(1-Q*Q)>-<0.5/Q><cos(Q)><Q*sqrt(1-Q*Q)_d>-<0.5/Q><sin(Q)_d><1/sqrt(1-Q*Q)>-<d_0.5><cos(Q)><sqrt(1-Q*Q)>+<0.5/Q><cos(Q)><d_Q*sqrt(1-Q*Q)>+<0.5/Q><d_sin(Q)><1/sqrt(1-Q*Q)>)";
    string ddyPolar="(<0.5_d><sin(Q)><sqrt(1-Q*Q)>-<0.5/Q><sin(Q)><Q*sqrt(1-Q*Q)_d>+<0.5/Q><cos(Q)_d><1/sqrt(1-Q*Q)>-<d_0.5><sin(Q)><sqrt(1-Q*Q)>+<0.5/Q><sin(Q)><d_Q*sqrt(1-Q*Q)>-<0.5/Q><d_cos(Q)><1/sqrt(1-Q*Q)>)";
    standardOperators["D/DX"]["Rn1.Phi1.Eta1.Rn2.Phi2.Eta2"]="(<1><1><1>"+ddxPolar+")$BAND.0.0.0.0.4.2 +("+ddxPolar+"<1><1><1>)$BAND.0.4.2.0.0.0";
    standardOperators["D/DY"]["Rn1.Phi1.Eta1.Rn2.Phi2.Eta2"]="(<1><1><1>"+ddyPolar+")$BAND.0.0.0.0.4.2 +("+ddyPolar+"<1><1><1>)$BAND.0.4.2.0.0.0";

    standardOperators["D/DX"]["Rn1.Phi1.Eta1"] = standardOperators["D/DX"]["Rn2.Phi2.Eta2"] =
            standardOperators["D/DX"]["Rn.Phi.Eta"]="(<0.5_d><cos(Q)><sqrt(1-Q*Q)>-<0.5/Q><cos(Q)><Q*sqrt(1-Q*Q)_d>-<0.5/Q><sin(Q)_d><1/sqrt(1-Q*Q)>"
            +string(                        "-<d_0.5><cos(Q)><sqrt(1-Q*Q)>+<0.5/Q><cos(Q)><d_Q*sqrt(1-Q*Q)>+<0.5/Q><d_sin(Q)><1/sqrt(1-Q*Q)>)$BAND.0.4.2");
    standardOperators["D/DY"]["Rn1.Phi1.Eta1"] = standardOperators["D/DY"]["Rn2.Phi2.Eta2"] =
            standardOperators["D/DY"]["Rn.Phi.Eta"]
            =       "(<0.5_d><sin(Q)><sqrt(1-Q*Q)>-<0.5/Q><sin(Q)><Q*sqrt(1-Q*Q)_d>+<0.5/Q><cos(Q)_d><1/sqrt(1-Q*Q)>"
            +string("-<d_0.5><sin(Q)><sqrt(1-Q*Q)>+<0.5/Q><sin(Q)><d_Q*sqrt(1-Q*Q)>-<0.5/Q><d_cos(Q)><1/sqrt(1-Q*Q)>)$BAND.0.4.2");
    standardOperators["D/DZ"]["Rn1.Phi1.Eta1"]=
            standardOperators["D/DZ"]["Rn2.Phi2.Eta2"] =
            standardOperators["D/DZ"]["Rn.Phi.Eta"]="(<0.5_d><1><Q>+<0.5/Q><1><(1-Q*Q)_d>"
            +string(                                "-<d_0.5><1><Q>-<0.5/Q><1><d_(1-Q*Q)>)$BAND.0.0.2");
    standardOperators["D/DX"]["Vec.Rn1.Phi1.Eta1"] =
            standardOperators["D/DX"]["Vec.Rn2.Phi2.Eta2"]
            =       "(<1><0.5_d><cos(Q)><sqrt(1-Q*Q)>-<1><0.5/Q><cos(Q)><Q*sqrt(1-Q*Q)_d>-<1><0.5/Q><sin(Q)_d><1/sqrt(1-Q*Q)>"
            +string("-<1><d_0.5><cos(Q)><sqrt(1-Q*Q)>+<1><0.5/Q><cos(Q)><d_Q*sqrt(1-Q*Q)>+<1><0.5/Q><d_sin(Q)><1/sqrt(1-Q*Q)>)$BAND.0.0.4.2");
    standardOperators["D/DY"]["Vec.Rn1.Phi1.Eta1"] =
            standardOperators["D/DY"]["Vec.Rn2.Phi2.Eta2"]
            =       "(<1><0.5_d><sin(Q)><sqrt(1-Q*Q)>-<1><0.5/Q><sin(Q)><Q*sqrt(1-Q*Q)_d>+<1><0.5/Q><cos(Q)_d><1/sqrt(1-Q*Q)>"
            +string("-<1><d_0.5><sin(Q)><sqrt(1-Q*Q)>+<1><0.5/Q><sin(Q)><d_Q*sqrt(1-Q*Q)>-<1><0.5/Q><d_cos(Q)><1/sqrt(1-Q*Q)>)$BAND.0.0.4.2");
    standardOperators["D/DZ"]["Vec.Rn1.Phi1.Eta1"] =
            standardOperators["D/DZ"]["Vec.Rn2.Phi2.Eta2"]="(<1><0.5_d><1><Q>+<1><0.5/Q><1><(1-Q*Q)_d>"
            +string(                       " -<1><d_0.5><1><Q>-<1><0.5/Q><1><d_(1-Q*Q)>)$BAND.0.0.0.2");

    standardOperators["D/DZ"]["Ion.Rn.Phi.Eta"]="(<1><0.5_d><1><Q>+<1><0.5/Q><1><(1-Q*Q)_d>"
            +string(                       " -<1><d_0.5><1><Q>-<1><0.5/Q><1><d_(1-Q*Q)>)$BAND.0.0.0.2+<D/Dz><1><1><1>")
            +string(                       "+<allOnes><NONORTHz_vg><allOnes><allOnes>");
    standardOperators["D/DZ"]["Vec.Ion"]="<1><D/Dz>";

    standardOperators["D/DZ"]["Rn1.Phi1.Eta1.Rn2.Phi2.Eta2"]="(<0.5_d><1><Q><1><1><1>+<0.5/Q><1><(1-Q*Q)_d><1><1><1>"
            +string(                       " -<d_0.5><1><Q><1><1><1>-<0.5/Q><1><d_(1-Q*Q)><1><1><1>)$BAND.0.0.2.0.0.0")
            +string(                           " +(<1><1><1><0.5_d><1><Q>+<1><1><1><0.5/Q><1><(1-Q*Q)_d>" )
            +string(                       " -<1><1><1><d_0.5><1><Q>-<1><1><1><0.5/Q><1><d_(1-Q*Q)>)$BAND.0.0.0.0.0.2");
    standardOperators["D/DZ"]["nSurface.Phi1.Eta1.Rn2.Phi2.Eta2"]="(<1><1><1><0.5_d><1><Q>+<1><1><1><0.5/Q><1><(1-Q*Q)_d>"
            +string(                       " -<1><1><1><d_0.5><1><Q>-<1><1><1><0.5/Q><1><d_(1-Q*Q)>)$BAND.0.0.0.0.0.2");
    standardOperators["D/DZ"]["Rn1.Phi1.Eta1.nSurface.Phi2.Eta2"]="(<0.5_d><1><Q><1><1><1>+<0.5/Q><1><(1-Q*Q)_d><1><1><1>"
            +string(                       " -<d_0.5><1><Q><1><1><1>-<0.5/Q><1><d_(1-Q*Q)><1><1><1>)$BAND.0.0.2.0.0.0");

    standardOperators["Z"]["Rn1.Phi1.Eta1.Rn2.Phi2.Eta2"]="(<Q><1><Q><1><1><1>+<1><1><1><Q><1><Q>)$BAND.0.0.2.0.0.2";

    // dipole, dipole velocity, dipole acceleration
    string potDer("<1/(Q*Q)>");
    if(smoothFunction!="1")potDer="<"+smoothFunction+"/(Q*Q)-"+smoothDerivative+"/Q>";
    standardOperators["DipAccX"]["Rn.Phi.Eta"]=potDer+"<cos(Q)><sqrt(1-Q*Q)>";
    standardOperators["DipAccY"]["Rn.Phi.Eta"]=potDer+"<sin(Q)><sqrt(1-Q*Q)>";
    standardOperators["DipAccZ"]["Rn.Phi.Eta"]=potDer+"<1><Q>";

    standardOperators["PotSolid"]["Z"]="<PotSolid>";
    standardOperators["PotSolid"]["Vec.Z"]="<1><PotSolid>";


    standardOperators["DipAccX"]["Vec.X"]=standardOperators["DipAccY"]["Vec.Y"]=standardOperators["DipAccZ"]["Vec.Z"]="<-1><Derivative_PotSolid>";
    standardOperators["DipAccX"]["X"]=standardOperators["DipAccY"]["Y"]=standardOperators["DipAccZ"]["Z"]="-<Derivative_PotSolid>";

    standardOperators["D/DX"]["X"]=standardOperators["D/DY"]["Y"]=standardOperators["D/DZ"]["Z"]="(<d_0.5>-<0.5_d>)";
    standardOperators["X"]["X"]=standardOperators["Y"]["Y"]=standardOperators["Z"]["Z"]="(<Q>)";


    // the mixed gauge dipole term
    string r0("chi[Rg,infty]*0.5");                   // factor 1/2 for anti-symmetrization included here
    string r1("chi[Rg,infty]*((Q-Rg)/(Q*Q))*0.5");    // factor 1/2 for anti-symmetrization included here
    string r2("chi[Rg,infty]*Rg*(2*Q-Rg)/(Q*Q)*0.5"); // factor 1/2 from kinetic energy
    string r3("chi[0,Rg]*Q+chi[Rg,infty]*Rg");
    string eta1("(1-Q*Q)");
    string eta2("sqrt"+eta1);

    string mixDip;
    // field-terms
    mixDip+= "LaserFz[t]<"+r3+"><1><Q>$BAND.0.0.2";
    mixDip+="+LaserFx[t]<"+r3+"><cos(Q)><"+eta2+">$BAND.0.4.2";
    mixDip+="+LaserFy[t]<"+r3+"><sin(Q)><"+eta2+">$BAND.0.4.2";

    // linear A-terms (anti-symmetrized)
    mixDip+="+iLaserAx[t](<"  +r0+"_d>< cos(Q)><"+eta2+">+<"+r1+">(<cos(Q)>< -Q*"+eta2+"_d>+<  sin(Q)_d><-1/"+eta2+">)"
            +           "+<d_"+r0+  "><-cos(Q)><"+eta2+">+<"+r1+">(<cos(Q)><d_Q*"+eta2+  ">+<d_sin(Q)  >< 1/"+eta2+">))$BAND.0.4.2";

    mixDip+="+iLaserAy[t](<"  +r0+"_d>< sin(Q)><"+eta2+">+<"+r1+">(<sin(Q)>< -Q*"+eta2+"_d>+<  cos(Q)_d>< 1/"+eta2+">)"
            +           "+<d_"+r0+  "><-sin(Q)><"+eta2+">+<"+r1+">(<sin(Q)><d_Q*"+eta2+  ">+<d_cos(Q)  ><-1/"+eta2+">))$BAND.0.4.2";

    mixDip+="+iLaserAz[t](<"       +r0+"_d><1><Q>+< "+r1+"><1><"       +eta1+"_d>"
            +           "+<d_(-1)*"+r0+"><1><Q>+<"   +r1+"><1><d_(-1)*"+eta1+">)$BAND.0.0.2";

    // quadratic A-terms...
    mixDip+="+LaserAsq[t]<-1+chi[Rg,infty]*pow[2]((Q-Rg)/Q)><0.5><1>$BAND.0.2.2"; // note: for Rg=0 no Asq-term

    mixDip+="+LaserAxx[t]<"+r2+"><cos(Q)*cos(Q)><"+eta1+">$BAND.0.4.4";
    mixDip+="+LaserAyy[t]<"+r2+"><sin(Q)*sin(Q)><"+eta1+">$BAND.0.4.4";
    mixDip+="+LaserAzz[t]<"+r2+"><1><Q*Q>$BAND.0.0.4";

    mixDip+="+LaserAxy[t]<"+r2+"><2*sin(Q)*cos(Q)><"+eta1+">$BAND.0.6.4";
    mixDip+="+LaserAxz[t]<"+r2+"><2*cos(Q)><Q*"+eta2+">$BAND.0.4.4";
    mixDip+="+LaserAyz[t]<"+r2+"><2*sin(Q)><Q*"+eta2+">$BAND.0.4.4";

    standardOperators["MixedGaugeDipole"]["Rn.Phi.Eta"]=unBracket(mixDip,"");

    // NOTE: a constraint for the first factor must be inserted
    standardOperators["MixedGaugeDipole"]["Ion.Rn.Phi.Eta"]=addConstraint("<1>("+mixDip+")","BAND",0,"0");
    standardOperators["MixedGaugeDipole"]["Ion.Rn.Phi.Eta"]+="+(LaserFx[t]<x>+LaserFy[t]<y>+LaserFz[t]<z>)<1><1><1>";
    standardOperators["MixedGaugeDipole"]["Ion.Rn.Phi.Eta"]+="+LaserFx[t]<allOnes><NONORTHx><allOnes><allOnes>";
    standardOperators["MixedGaugeDipole"]["Ion.Rn.Phi.Eta"]+="+LaserFy[t]<allOnes><NONORTHy><allOnes><allOnes>";
    standardOperators["MixedGaugeDipole"]["Ion.Rn.Phi.Eta"]+="+LaserFz[t]<allOnes><NONORTHz><allOnes><allOnes>";

    // quadratic antisym term
    standardOperators["MixedGaugeDipole"]["Ion.Rn.Phi.Eta"]+="-LaserAsqHalf[t]<allOnes><NONORTH1><allOnes><allOnes>";

    // pure length gauge coupling of ionic and states
    standardOperators["MixedGaugeDipole"]["Neutral"]="(LaserFx[t]<x>+LaserFy[t]<y>+LaserFz[t]<z>)";
    standardOperators["MixedGaugeDipole"]["Vec.Ion"]="<1>(LaserFx[t]<x>+LaserFy[t]<y>+LaserFz[t]<z>)";

    // remove Asq terms from Neutral
    standardOperators["MixedGaugeDipole"]["Neutral"]+="-LaserAsqHalf[t]<1>";
    //standardOperators["MixedGaugeDipole"]["Vec.Ion"]="<1>(-<1>)";

    // commutators - compose of basic 1d and polar
    standardOperators        ["Commutator"]["ValDerX"]=
            standardOperators["Commutator"]["ValDerY"]=
            standardOperators["Commutator"]["ValDerZ"]=
            standardOperators["Commutator"]["ValDerX1"]=
            standardOperators["Commutator"]["ValDerX2"]=
            "-1/2<0,1>+1/2<1,0>+iLaserAz[t]<0,0>";

    standardOperators         ["Commutator"]["ValDerX.Y"]
            =standardOperators["Commutator"]["ValDerX.kY"]
            =standardOperators["Commutator"]["ValDerY.X"]
            =standardOperators["Commutator"]["ValDerY.kX"]
            ="("+standardOperators["Commutator"]["ValDerX"]+")<Id>";

    standardOperators         ["Commutator"]["ValDerX1.X2"]
            =standardOperators["Commutator"]["ValDerX1.kX2"]
            =standardOperators["Commutator"]["ValDerX2.X1"]
            =standardOperators["Commutator"]["ValDerX2.kX1"]
            ="("+standardOperators["Commutator"]["ValDerX1"]+")<Id>";

    standardOperators         ["Commutator"]["ValDerX.Z"]
            =standardOperators["Commutator"]["ValDerX.kZ"]
            =standardOperators["Commutator"]["ValDerZ.X"]
            =standardOperators["Commutator"]["ValDerZ.kX"]
            ="("+standardOperators["Commutator"]["ValDerX"]+")<Id>";

    standardOperators         ["Commutator"]["ValDerX.Y.Z"]
            =standardOperators["Commutator"]["ValDerX.kY.Z"]
            =standardOperators["Commutator"]["ValDerX.Y.kZ"]
            =standardOperators["Commutator"]["ValDerX.kY.kZ"]
            =standardOperators["Commutator"]["ValDerX.Z.Y"]
            =standardOperators["Commutator"]["ValDerX.kZ.Y"]
            =standardOperators["Commutator"]["ValDerX.Z.kY"]
            =standardOperators["Commutator"]["ValDerX.kZ.kY"]
            =standardOperators["Commutator"]["ValDerY.X.Z"]
            =standardOperators["Commutator"]["ValDerY.kX.Z"]
            =standardOperators["Commutator"]["ValDerY.X.kZ"]
            =standardOperators["Commutator"]["ValDerY.kX.kZ"]
            =standardOperators["Commutator"]["ValDerY.Z.X"]
            =standardOperators["Commutator"]["ValDerY.kZ.X"]
            =standardOperators["Commutator"]["ValDerY.Z.kX"]
            =standardOperators["Commutator"]["ValDerY.kZ.kX"]
            =standardOperators["Commutator"]["ValDerZ.X.Y"]
            =standardOperators["Commutator"]["ValDerZ.kX.Y"]
            =standardOperators["Commutator"]["ValDerZ.X.kY"]
            =standardOperators["Commutator"]["ValDerZ.kX.kY"]
            =standardOperators["Commutator"]["ValDerZ.Y.X"]
            =standardOperators["Commutator"]["ValDerZ.kY.X"]
            =standardOperators["Commutator"]["ValDerZ.Y.kX"]
            =standardOperators["Commutator"]["ValDerZ.kY.kX"]
            ="("+standardOperators["Commutator"]["ValDerX"]+")<Id><Id>";

    standardOperators["Commutator"]["Phi.Eta.ValDerRn"]=
            "iLaserAz[t]<1><Q><0,0>+iLaserAx[t]<cos(Q)><sqrt(1-Q*Q)><0,0>+iLaserAy[t]<sin(Q)><sqrt(1-Q*Q)><0,0>+1/2<1><1><1,0>-1/2<1><1><0,1>";

    standardOperators["Commutator"]["Phi1.Eta1.Phi2.Eta2.ValDerRn1.Rn2"]="iLaserAz[t]<1><Q><1><1><0,0><1>+iLaserAx[t]<cos(Q)><sqrt(1-Q*Q)><1><1><0,0><1>+iLaserAy[t]<sin(Q)><sqrt(1-Q*Q)><1><1><0,0><1>+1/2<1><1><1><1><1,0><1>-1/2<1><1><1><1><0,1><1>";
    standardOperators["Commutator"]["Phi1.Eta1.Phi2.Eta2.Rn1.ValDerRn2"]="<1><1>(iLaserAz[t]<1><Q><1><0,0>+iLaserAx[t]<cos(Q)><sqrt(1-Q*Q)><1><0,0>+iLaserAy[t]<sin(Q)><sqrt(1-Q*Q)><1><0,0>+1/2<1><1><1><1,0>-1/2<1><1><1><0,1>)";

    standardOperators["Commutator"]["Phi.Eta.ValDerRn"]=
            "<1><Q>iLaserAz[t]<0,0>+<cos(Q)><sqrt(1-Q*Q)>iLaserAx[t]<0,0>+<sin(Q)><sqrt(1-Q*Q)>iLaserAy[t]<0,0>+<1><1>1/2<1,0>-<1><1>1/2<0,1>";

    standardOperators         ["Commutator"]["Phi.Eta.ValDerRn.Vec"]
            =standardOperators["Commutator"]["Phi1.Eta1.ValDerRn1.Vec"]
            =standardOperators["Commutator"]["Phi2.Eta2.ValDerRn2.Vec"]
            ="("+standardOperators["Commutator"]["Phi.Eta.ValDerRn"]+")<Id>";

    standardOperators         ["Commutator"]["Phi1.Eta1.ValDerRn1.Phi2.Eta2.Rn2"]
            =standardOperators["Commutator"]["Phi1.Eta1.ValDerRn1.Phi2.Eta2.kRn2"]
            =standardOperators["Commutator"]["Phi2.Eta2.ValDerRn2.Phi1.Eta1.Rn1"]
            =standardOperators["Commutator"]["Phi2.Eta2.ValDerRn2.Phi1.Eta1.kRn1"]
            ="("+standardOperators["Commutator"]["Phi.Eta.ValDerRn"]+")<Id><Id><Id>";

    standardOperators["Commutator"]["Phi1.Eta1.Phi2.Eta2.ValDerRn1.Rn2"]="iLaserAz[t]<1><Q><1><1><0,0><1>+iLaserAx[t]<cos(Q)><sqrt(1-Q*Q)><1><1><0,0><1>+iLaserAy[t]<sin(Q)><sqrt(1-Q*Q)><1><1><0,0><1>+1/2<1><1><1><1><1,0><1>-1/2<1><1><1><1><0,1><1>";
    standardOperators["Commutator"]["Phi1.Eta1.Phi2.Eta2.Rn1.ValDerRn2"]="iLaserAz[t]<1><1><1><Q><1><0,0>+iLaserAx[t]<1><1><cos(Q)><sqrt(1-Q*Q)><1><0,0>+iLaserAy[t]<1><1><sin(Q)><sqrt(1-Q*Q)><1><0,0>+1/2<1><1><1><1><1><1,0>-1/2<1><1><1><1><1><0,1>";
    standardOperators["Commutator"]["Phi1.Eta1.Phi2.Eta2.kRn1.ValDerRn2"]="iLaserAz[t]<1><1><1><Q><1><0,0>+iLaserAx[t]<1><1><cos(Q)><sqrt(1-Q*Q)><1><0,0>+iLaserAy[t]<1><1><sin(Q)><sqrt(1-Q*Q)><1><0,0>+1/2<1><1><1><1><1><1,0>-1/2<1><1><1><1><1><0,1>";;
    standardOperators["Commutator"]["Phi1.Eta1.Phi2.Eta2.ValDerRn1.kRn2"]="iLaserAz[t]<1><Q><1><1><0,0><1>+iLaserAx[t]<cos(Q)><sqrt(1-Q*Q)><1><1><0,0><1>+iLaserAy[t]<sin(Q)><sqrt(1-Q*Q)><1><1><0,0><1>+1/2<1><1><1><1><1,0><1>-1/2<1><1><1><1><0,1><1>";

}


//string OperatorDefinition::undefined="UNDEFINED";

void OperatorDefinition::setParameters(const string Definition){
    size_t pos=Definition.find("<<",0);
    while(pos!=string::npos){
        extractParameters(Definition.substr(pos,Definition.find(">>",pos)+2-pos));
        pos=Definition.find("<<",pos+2);
    }
}

void OperatorDefinition::truncationRadius(double Rsmooth, double Rtrunc){
    Algebra::addSpecialConstant("Rsmooth",Rsmooth);
    Algebra::addSpecialConstant("Rtrunc",Rtrunc);
}

string OperatorDefinition::extractParameters(const string Name)
{
    if(tools::subStringCount(Name,"<<")>1)ABORT("only for single standard operator in definition");
    string name=tools::stringInBetween(Name,"<<",">>",true);
    size_t pos=name.find(":");
    if(pos==string::npos)return Name;

    vector<string>pars=tools::splitString(name.substr(pos+1),',');
    for(vector<string>::iterator par=pars.begin();par!=pars.end();par++){
        size_t sep=par->find("=");
        if(sep==string::npos)ABORT("invalid parameter specification in "+Name+"\n specify as in Name(par=123)");
        Algebra::addSpecialConstant(par->substr(0,sep),tools::string_to_double(par->substr(sep+1)));
    }

    return name.substr(0,pos);
}

std::string OperatorDefinition::get(std::string Name, std::string Hierarchy){
    string coor=Index::coordinates(Hierarchy);
    Name=extractParameters(Name); // e.g. for <<MixedGaugeDipole:Rg=xx>>: add value Rg to parameter list with value xx

    ReadInput::main.obsolete("Source","undoundDOF","do not use, instead specify region as Rn1, Rn2, Rn1.Rn2, etc.");

    if(not tools::hasKey(standardOperators,Name)){
        return "undefined <<"+Name
                +">>\n available operators:\n"+tools::listMapKeys(standardOperators,"\n");
    }

    std::vector<string> keys=tools::vectorMapKeys(standardOperators[Name]);
    std::vector<unsigned int> perm;
    string opDef="";
    for(unsigned int k=0;k<keys.size();k++){
        perm=tools::permutation(tools::splitString(keys[k],'.'),tools::splitString(coor,'.'));
        if(perm.size()!=0){
            opDef=unBracket(standardOperators[Name][keys[k]]);
            break;
        }
    }
    
    if(opDef.length()==0 and (Name=="Overlap" or Name=="1")){  // use default Overlap: <1><1>...
        std::vector<std::string>cov=tools::splitString(coor,'.');
        for(auto c: cov){
            // Hybrids may have overlap between blocks, except Subspace&Complement, which are orthogonal by construction
            if(c.find("&")!=std::string::npos and Index::getAxisSubset(c)!="Subspace&Complement")
                opDef+="<allOnes>";
            else if(c=="Ndim"){
                PrintOutput::DEVwarning("hack: assuming 3d for NDim");
                opDef+="<1><1><1>";
            }
            else
                opDef+="<1>";
        }
        standardOperators["Overlap"][coor]=opDef;
        opDef=unBracket(standardOperators["Overlap"][coor]);
    }

    if(opDef.length()==0 and (Name=="GridWeights")){  // use default Overlap: <1><1>...
        standardOperators[Name][coor]="<GridWeight>";
        for(int k=0;k<std::count(coor.begin(),coor.end(),'.');k++){
            standardOperators[Name][coor]+="<GridWeight>";
        }
        opDef=unBracket(standardOperators[Name][coor]);
    }

    if(opDef.length()==0)return "undefined <<"+Name+">> for "+coor
            +"\n available for permutations of\n"+tools::listMapKeys(standardOperators[Name],"\n");

    bool permuted=false;
    for(unsigned int i=0;i<perm.size();i++) if(perm[i]!=i) {permuted=true; break;}
    if(not permuted) return opDef;

    // for now, we cannot handle brackets outside factors
    if(tools::findFirstOutsideBrackets(opDef,"(","<",">")!=string::npos)
        ABORT("cannot handle brackets ouside factors\nDef: "+opDef);

    return opPermuted(opDef,perm);
}

string OperatorDefinition::addConstraint(string Def, string ConstType, int Pos, string C){
    size_t cpos=Def.find(ConstType);
    while(cpos!=string::npos){
        for(int k=0;k<Pos+1 and cpos!=string::npos;k++)cpos=Def.find(".",cpos);
        if(cpos==string::npos)
            ABORT(Str("Def="+Def+": no position=")+Pos+("in constraint $"+ConstType));
        Def.insert(cpos,"."+C);
        cpos=Def.find(ConstType,cpos);
    }
    return Def;
}

string OperatorDefinition::opPermuted(const string & Op, const std::vector<unsigned int> Perm){
    string oPerm;
    string op(Op);
    if(tools::subStringCount(op,"<")%Perm.size()!=0)
        ABORT("incorrect number of factors in operator definition\n"
              +Op+": factors="+tools::str(tools::subStringCount(op,"<"))+", dimension="+tools::str(Perm.size()));

    vector<string> term,sign;
    tools::splitString(Op,"+-",term,sign,"<|",">>");
    for(unsigned int k=0;k<term.size();k++){
        string front=sign[k]+term[k].substr(0,term[k].find('<'));
        string back=term[k].substr(term[k].find('<'));
        vector<string>fac=tools::splitString(back,'>');
        while(fac.back()=="")fac.pop_back(); // remove trailing empty strings
        oPerm+=front;
        // permute into oPerm
        for(unsigned int l=0;l<Perm.size();l++){
            oPerm+=fac[Perm[l]]+">";
        }
        // permute the constraints
        if(fac.size()>Perm.size()){
            vector<string>lCons=tools::splitString(fac.back(),'.');
            if(lCons.size()<Perm.size()+1)
                ABORT(Str("number of constraints in ")+term[k]+"smaller than permutation "+Perm);
            if(tools::cropString(lCons[0])!="$BAND")ABORT("only band qualifiers implemented, is: '"+lCons[0]+"'\n"+Op);
            oPerm+=tools::cropString(lCons[0]);
            for(unsigned int l=0;l<Perm.size();l++)oPerm+="."+lCons[1+Perm[l]];
        }
    }

    return oPerm;
}

void OperatorDefinition::syntax() const {
    if(find("<d_-")!=string::npos or find("-_d>")!=string::npos or find("<d_+")!=string::npos or find("+_d>")!=string::npos){
        //note: tests is too crude...
        ABORT("malformed operator: "+*this+"\ndo not use signs next to _d> or <d_, enclose in brackets, e.g. <d_(-0.5)> etc.");
    }
}

void OperatorDefinition::Test(){

    OperatorDefinition d;
    cout<<"Operators and some permutations"<<endl;
    cout<<"Laplacian(Rn.Phi.Eta): "<<get("Laplacian","Rn.Phi.Eta")<<endl;
    cout<<"Laplacian(Phi.Rn.Eta): "<<get("Laplacian","Phi.Rn.Eta")<<endl;
    cout<<"Laplacian(Phi.Rn.bad): "<<get("Laplacian","Phi.Eta.bad")<<endl;
    cout<<"Laplacian(Phi.Eta.Rn): "<<get("Laplacian","Phi.Eta.Rn")<<endl;
    cout<<"Laplacian(Rn1.Rn2.Phi1.Phi2.Eta1.Eta2): "<<get("Laplacian","Rn1.Rn2.Phi1.Phi2.Eta1.Eta2")<<endl;
    cout<<"D/DX(Rn.Phi.Eta): "<<get("D/DX","Rn.Phi.Eta")<<endl;
    cout<<"D/DX(Phi.Rn.Eta): "<<get("D/DX","Phi.Rn.Eta")<<endl;
    cout<<"D/DX(Eta.Phi.Rn): "<<get("D/DX","Eta.Phi.Rn")<<endl;
}
OperatorDefinition OperatorDefinition::first(bool Unsign) const{

    OperatorDefinition res;
    //    if(not OperatorData::isStandard(*this))return *this;
    if(not isStandard(0,0))return *this;
    res=*this;
    if(Unsign){
        vector<string>elem,seps;
        tools::splitString(res,"+-",elem,seps,"(<",")>");
        res=elem[0];
    }
    // first factor must be a single-bracket expression, else return unchanged
    return res.find(">")<res.find("<<") and res.find(">")<res.find("[[")?OperatorDefinition(res.substr(0,res.find(">")+1)):*this;
}

std::vector<OperatorDefinition> OperatorDefinition::factors() const{
    OperatorDefinition fac;
    OperatorDefinition remain(*this);
    std::vector<OperatorDefinition>res;
    while(""!=(fac=remain.firstFactor())){
        res.push_back(remain.substr(0,fac.length()));
        remain=remain.substr(fac.length());
    }
    return res;
}

std::vector<OperatorDefinition> OperatorDefinition::terms() const {
    OperatorDefinition udef=unBracket(*this);
    vector<string> term,sign,ret0;
    std::vector<OperatorDefinition> ret;
    OperatorDefinition def(tools::cropString(udef));

    if(def.substr(0,1)=="'"){
        if(def.find("'",1)!=def.length()-1)ABORT("unmatched first ' in "+udef);
        def=def.substr(1);
        def=def.substr(0,def.length()-1);
    }

    // strip overall round brackets
    if(def[0]=='(' and def[def.length()-1]==')')def=def.substr(1,def.length()-2);

    tools::splitString(udef,"+-",term,sign,"<|(",">>)");
    for(unsigned int n=0;n<term.size();n++)
        if(term[n].length()!=0 and term[n]!=" ")ret0.push_back(tools::cropString(sign[n]+term[n]));

    // un-bracket grouped terms
    for(unsigned int k=0;k<ret0.size();k++){
        if(ret0[k].rfind(")")==ret0[k].length()-1){
            string par=tools::cropString(ret0[k].substr(0,ret0[k].find("(")));
            tools::splitString(tools::stringInBetween(ret0[k],"(",")"),"+-",term,sign,"<|(",">>)");

            for(unsigned int l=0;l<term.size();l++){
                // adjust the sign
                string spar=par;
                switch (spar[0]){
                case '+':
                    if(sign[l]=="-")spar[0]='-';
                    break;
                case '-':
                    if(sign[l]=="-")spar[0]='+';
                    break;
                default:
                    spar=sign[l]+par;
                }
                ret.push_back(spar+term[l]);
            }
        } else {
            ret.push_back(ret0[k]);
        }
    }

    for(unsigned int k=0;k<ret.size();k++)
        if(ret[k].find(")")==ret[k].length()-1)
            ABORT("for now, cannot have nested (...) in operator: "+udef);

    return ret;
}

OperatorDefinition OperatorDefinition::extractBlock(size_t I, size_t J) const {
    OperatorDefinition extract;
    //    vector<string> trms(OperatorData::terms(*this));
    vector<OperatorDefinition> trms(terms());
    for(size_t k=0;k<trms.size();k++)
    {

        // adjust constraints (if any)
        if(trms[k].find("$")!=string::npos){
            // first constraint parameter and remove
            size_t pos0=trms[k].find(".",trms[k].find("$"));
            size_t pos1=min(trms[k].length(),trms[k].find(".",pos0+1));
            trms[k].erase(pos0,pos1-pos0);
        }

        // loop through terms
        //        string blk=OperatorData::first(trms[k]); // first factor of term
        string blk=trms[k].first(); // first factor of term
        size_t posL=blk.find("<");
        size_t posR=blk.find(">");
        size_t posM=blk.find(",");
        if(posM==string::npos and blk.find("<1>")==string::npos and  blk.find("<allOnes>")==string::npos)
            ABORT("term["+tools::str(k)+"] = "+trms[k]+" of "+*this+" is not block-index: "
                  +blk+" must of the form <i,j>, <1>, or <allOnes>");

        // if matches <I,J> (or <1> or <allOnes>), remove first factor and append to block definition
        if(blk.find("<allOnes>")==posL or (I==J and blk.find("<1>")==posL))
            extract+=trms[k].substr(0,posL)+trms[k].substr(posR+1);
        else {
            auto si = blk.substr(posL+1,posM-posL-1);
            auto sj = blk.substr(posM+1,posR-posM-1);
            auto ci = tools::string_to_int(si);
            auto cj = tools::string_to_int(sj);
            if(posM != string::npos && ci ==int(I) and cj==int(J)){
                extract+=trms[k].substr(0,posL)+trms[k].substr(posR+1);
            }
        }

    }
    return extract;
}

static std::string originalHierarchy(const std::string Hierarchy){
    // replace k- and spec-prefixes in Hierarchy
    std::vector<std::string> coors=tools::splitString(Hierarchy,'.');
    for(auto &c: coors){
        if(c[0]=='k')c=c.substr(1);
        if(c.substr(0,4)=="spec")c=c.substr(4);
    }
    return tools::joinString(coors,".");
}

std::vector<OperatorDefinition> OperatorDefinition::singleTerms(const std::string & IHierarchy, const std::string & JHierarchy) const {

    OperatorDefinition def=OperatorDefinition::unBracket(*this);

    // expand if on equal coordinates and not hybrid
    if(Index::coordinates(IHierarchy)==Index::coordinates(JHierarchy) and IHierarchy.find("&")==std::string::npos)
        def=def.expandStandard(originalHierarchy(IHierarchy));

    // new brackets may have appeared - unbracket again
    return unBracket(def).terms();
}

// expand a single predefined operator that is embedded within single-factor terms
OperatorDefinition OperatorDefinition::expandStandard(std::string Hierarchy) const {
    // nothing to be expanded
    if(find("<<")==string::npos)return *this;

    // create the table of standard operators
    //    OperatorData::setStandard();
    OperatorDefinition::setup();

    vector<string> term,sep;
    tools::splitString(OperatorDefinition::unBracket(*this),"+-",term,sep,"<(",">)");

    vector<string> axis=tools::splitString(Hierarchy,'.');

    // determine floor axis kax - either second appearance or non-axis
    int kax=0;
    for(;kax<int(axis.size());kax++){
        int l=0;
        for(;l<kax;l++){
            if(axis[l]==axis[kax])break;
        }
        if(l!=kax or axis[kax]=="NONE")break;
    }
    axis.resize(kax);

    string def;
    for(unsigned int t=0;t<term.size();t++){
        if(tools::subStringCount(term[t],"<<")==0){
            def+=sep[t]+term[t];
        }
        else if(tools::subStringCount(term[t],"<<")==1){
            unsigned int kstart=tools::subStringCount(term[t].substr(0,term[t].find("<<")),">");
            unsigned int kend=axis.size()-tools::subStringCount(term[t].substr(term[t].find(">>")),"<");
            string coor=axis[kstart];
            for(unsigned int k=kstart+1;k<kend;k++)coor+="."+axis[k];
            string op=OperatorDefinition::extractParameters(tools::stringInBetween(term[t],"<<",">>",true));
            if(op=="")ABORT("empty standard operator in "+term[t]);

            string expand=OperatorDefinition::get(op,coor);
            // embed possible constraints
            size_t icon=expand.find("$");
            while(icon!=string::npos){
                for(unsigned int k=0;k<kstart;k++)expand.insert(expand.find(".",icon),".*");
                if(kend!=axis.size())ABORT("for now, not constraint extension to lower levels");
                icon=expand.find("$",icon+1);
            }

            def+=sep[t]+term[t].substr(0,term[t].find("<<"))+"("+expand+")";
        }
        else ABORT("cannot multiply standard operators: "+term[t]);

    }
    if(ReadInput::main.flag("DEBUGshowOp","show the expanded operator definitions"))PrintOutput::lineItem(*this+" expands to: ",def);
    return def;
}


using namespace std;

static std::map<std::string,std::string> predefinedInputs={
    {"<<LaserVelocityZX>>","iLaserAz[t]<<D/DZ>>+iLaserAx[t]<<D/DX>>"}
    ,{"<<LaserVelocityRotatingFrame>>","iLaserAabs[t]<<D/DZ>>+iLaserAdTheta[t]<<CoriolisTheta>>"}
    ,{"<<LaserLengthZX>>","LaserFz[t]<<Z>>+LaserFx[t]<<X>>-0.5*LaserAsq[t]<<Overlap>>"}
    ,{"<<LaserLengthRotatingFrame>>","LaserFabs[t]<<Z>>-0.5*LaserAsq[t]<<Overlap>>+iLaserFdTheta[t]<<CoriolisTheta>>"}
    ,{"<<DipVelX>>","i<<D/DX>>"}
    ,{"<<DipVelY>>","i<<D/DY>>"}
    ,{"<<DipVelZ>>","i<<D/DZ>>"}
    ,{"<<DipLenX>>","<<X>>"}
    ,{"<<DipLenY>>","<<Y>>"}
    ,{"<<DipLenZ>>","<<Z>>"}
};

static bool isPredefined(std::string Def){
    if(Def=="<<Overlap>>")return true;
    Def=tools::stringInBetween(Def,"<<",">>",true);
    for(auto &d: predefinedInputs)
        if(d.first=="<<"+Def+">>")return true;
    for(auto &d: OperatorDefinition::standardOperators)
        if(d.first==OperatorDefinition::extractParameters(Def))return true;
    return false;

}
static std::string availableDefinitions(){
    std::string res;
    for(auto &d: predefinedInputs){
        res+=d.first+": "+d.second+"\n";
    }
    for(auto &d: OperatorDefinition::standardOperators){
        res+=d.first+" for "+tools::listMapKeys(d.second)+"\n";
    }
    res.pop_back();
    return res;
}

OperatorDefinition::OperatorDefinition(const std::string InputDef, std::string Hierarchy, bool Abort)
    :OperatorDefinition(InputDef)
{
    _abort=Abort;
    if(InputDef=="")return;
    if(Hierarchy.find("&")!=std::string::npos)return; // do not expand operators on hybrid axis

    std::string Def=predefinedInputs.count(InputDef)?predefinedInputs[InputDef]:InputDef;

    if(Def.find("<<1>>")!=string::npos)
        ABORT("OBSOLETE operator definition <<1>>, replace by <<Overlap>>");

    if(tools::subStringCount(Def,">>")!=tools::subStringCount(Def,"<<"))
        ABORT("unmatched <<...>> in "+Def);

    // special operators are not resolved
    if(Def.find("[[")<Def.find("<")){
        assign(Def);
        return;
    }

    // remove hybrid axes
    std::string hierarchy;
    hierarchy=Hierarchy;

    // parent Hierarchy and original input coordinates
    string parentHier=parentHierarchy(hierarchy);
    string inputCoor=inputCoordinates(parentHier);


    // split into terms, apply to each term
    string def;
    std::vector<OperatorDefinition>term=OperatorDefinition(Def).terms();

    // commutator for hybrid problem
    if(Def=="<<Commutator>>" and inputCoor.find('&')<inputCoor.find('.'))term[0]=std::string("<1><<Commutator>>");


    // empty operator - do nothing
    if(term[0]=="");

    // multiple terms
    else if(term.size()>1){
        for(size_t k=0;k<term.size();k++)
            def+=OperatorDefinition(term[k],inputCoor);
    }

    // factor<singleFactor>OpStrX
    else if(""!=firstSingleFactor(term[0])){
        string parF=parameter(term[0],"+");
        string singF=firstSingleFactor(term[0]);
        string debug=OperatorDefinition(remainder(term[0],singF),inputCoor.substr(inputCoor.find(".")+1));
        std::vector<OperatorDefinition>termY=OperatorDefinition(debug).terms();
        for(size_t k=0;k<termY.size();k++)
            def+=parameter(termY[k],parF)+compose(noParameter(singF),noParameter(termY[k]));
        // ugly hack
        if(OperatorDefinition(debug).terms().size()==0)
            def=parameter("",parF)+compose(noParameter(singF),noParameter(""));

    }

    // factor<<PredefinedOp>>OpStrX
    else if(""!=firstPredefined(term[0])){
        if(term[0].find("<<")>term[0].find(">>"))
            ABORT("for now, cannot have products of pre-defined operators, is: "+Def+" of "+inputCoor);

        // get the remainder
        string predef=firstPredefined(term[0]);
        if(predef.find("unresolved")!=std::string::npos){assign(def);return;}
        std::vector<OperatorDefinition>termX=OperatorDefinition(remainder(term[0],predef)).terms();
        if(not termX.size())termX.push_back(OperatorDefinition());// a really ugly historical hack

        // remove trailing coordinates for expanding predef
        string predefCoor=inputCoor;
        for(size_t posF=termX[0].find("<");posF!=string::npos;posF=termX[0].find("<",posF+1)){
            predefCoor=predefCoor.substr(0,predefCoor.rfind('.')-1);
        }
        if(predefCoor=="")ABORT("too few coordinates for expanding "+Def+" given "+inputCoor);

        string param=parameter(term[0],"+");
        vector<OperatorDefinition>termP=OperatorDefinition().construct(predef.substr(2,predef.length()-4),predefCoor).terms();
        std::string res=termP.size()?termP[0]:OperatorDefinition();
        while(res.find_first_of("+- ")==0)res=res.substr(1);
        if(predef.find(res)==2){
            def="unresolved: "+predef;
            if(_abort)ABORT(availableDefinitions()
                                     +"\navailable predefined operators above,"
                                     +"\nfailed to resolve "+predef+" for "+predefCoor+" derived from "+inputCoor);
        }
        for(size_t k=0;k<termP.size();k++){
            string parP=parameter(termP[k],param);
            for(size_t l=0;l<termX.size();l++){
                def+=parameter(termX[l],parP)+compose(noParameter(termP[k]),noParameter(termX[l]));
            }
        }
    }

    // OpStrX
    else {
        def=OperatorDefinition(term[0],inputCoor);
    }

    // drop terms that have factors !=<1> on certain coordinates
    //if(Drop)def=OperatorDefinitionNew(def,inputCoor).dropTerms(hierarchy);

    // on final level, permute terms for floor on bottom
    assign(def);

    if(inputCoor!=Index::coordinates(parentHier))
        assign(permuteFEM(parentHier));
}

std::string OperatorDefinition::defaultOverlap(const Index *Idx){
    std::string gridOv;
    for (const Index *ix = Idx; ix != 0; ix = ix->descend())
        if (not OperatorDefinition::skipAxis(ix))
            gridOv += ix->basis()->grid() ? "<GridWeight>" : "<1>";
    return gridOv;
}

std::string OperatorDefinition::noSign(std::string Term){
    if(Term.find_first_not_of(" ")==string::npos)return "";
    Term=Term.substr(Term.find_first_not_of(" "));
    if(Term.find_first_of("+-")!=0)return Term;
    return Term.substr(1);
}

/// extract parameter para from Term, if Param!="" return product with proper sign
std::string OperatorDefinition::parameter(std::string Term, std::string Param){
    string par=Term.substr(0,Term.find_first_of("<("));
    if(par==Term)par=Term.substr(0,Term.find("[["));

    if (par.find("[t]")!=string::npos){
        Algebra* alg=new Algebra(par.substr(int(par[0]=='+' or par[0]=='-')));
        if(alg->isAlgebra())Parameters::add(alg);
    }

    string sig=sign(par,sign(Param,"+"));
    string par0=noSign(par),Param0=noSign(Param);
    if(par0=="" or Param0=="")
        return   sig+Param0+par0;
    else
        return par="("+sig+Param0+"*"+par0+")";
}

/// remove sign from Par and return sign (no sign is returns +)
std::string OperatorDefinition::sign(std::string Par, const std::string OtherSign){
    if(OtherSign!="+" and OtherSign!="-")ABORT("OtherSign must be + or -, is \""+OtherSign+"\"");

    if(Par.find_first_not_of(" ")==string::npos)return OtherSign;
    Par=Par.substr(Par.find_first_not_of(" "));
    string s="+";
    if(Par.find_first_of("+-")==0)s=Par.substr(0,1);
    if(OtherSign==s)return "+";
    return "-";
}

std::string OperatorDefinition::remainder(std::string Term,std::string First){
    // remainder starts right after First (ignore $-constraint in First)
    string first=First.substr(0,First.find("$"));
    string rem=Term.substr(Term.find(first)+first.length());
    if(rem.find(">")<rem.find("$"))return rem; // actual remainder, not just constraint
    if(rem.find("]]")<rem.find("$"))return rem; // in case remainder is "special op"
    return string();
}

std::string OperatorDefinition::noParameter(std::string Term){
    // noParameter starts at first of (,<,$,[[
    size_t r=min(Term.find_first_of("(<$"),Term.find("[["));
    if(r==string::npos)return ""; // no remainder
    return Term.substr(r);
}

OperatorDefinition OperatorDefinition::firstFactor() const {

    OperatorDefinition res;
    if(""!=(res=firstSingleFactor(*this)))return res;
    if(""!=(res=firstPredefined(*this)))return res;
    if(""!=(res=firstSpecial()))return res;
    if(""!=(res=firstBracketed()))return res;
    return std::string("");
}

std::string OperatorDefinition::firstSingleFactor(const std::string Term) const {
    // the first factor is the first "<factorDev>" string
    // that is not preceded by >,]], or unclosed (

    string fac=Term.substr(0,Term.find(">",Term.find_first_of("<"))+1);
    if(fac=="")return ""; // no factor at all
    if(fac.find("<<")!=string::npos)return "";      // first is <<predefined>>
    if(fac.find("<")==string::npos)ABORT("bad factor "+fac+" in "+Term);
    if(     std::count(fac.begin(),fac.begin()+fac.find("<"),'(')!=
            std::count(fac.begin(),fac.begin()+fac.find("<"),')'))return ""; // is inside brackets
    if(fac.find("]]")!=string::npos)return "";  // preceeded by special operator

    // remove possible parameters
    fac=fac.substr(fac.find("<"));

    // include the first constraint (if any)
    if(fac!="" and Term.find("$")!=string::npos){
        string con=Term.substr(Term.find("$"));
        fac=fac+con.substr(0,con.find(".",con.find(".")+1));
    }
    return fac;
}

std::string OperatorDefinition::firstPredefined(std::string Term) const {
    // the first factor is the first "<<predefined>>" string
    // that is not preceded by >,]], or unclosed (

    if(Term.find("<")==string::npos){
        if(Term.find("[[")!=string::npos)return Term;
        if(_abort)ABORT("cannot find operator in "+Term);
        return "unresolved: "+Term;
    }
    if(Term.find("<<")==string::npos)return "";
    string fac=Term.substr(0,Term.find("<<"));
    if(     std::count(fac.begin(),fac.end(),'(')!=
            std::count(fac.begin(),fac.end(),')'))return ""; // is inside brackets
    if(fac.find("]]")!=string::npos or fac.find(">")!=string::npos)return "";
    fac=Term.substr(Term.find("<<"));
    fac=fac.substr(0,Term.find(">>")+2);
    if(not isPredefined(fac)){
        if(not _abort)return "unresolved: "+fac;
        ABORT("\nAvailable predefined operators:\n\n"+
                                   availableDefinitions()+"\n\n"+fac
                                   +" does not exist, available definitions above");
    }

    return fac;
}

std::string OperatorDefinition::firstSpecial() const {
    // the first factor is the first "<<predefined>>" string
    // that is not preceded by >,]], or unclosed (

    auto beg=find("[[");
    if(beg>find_first_of(">"))return "";
    std::string fac=substr(beg,find("]]")+2);
    if(     std::count(fac.begin(),fac.end(),'(')!=
            std::count(fac.begin(),fac.end(),')'))return ""; // is inside brackets
    return fac;
}
std::string OperatorDefinition::firstBracketed() const {
    auto beg=std::min(find('<'),find("[["));
    int open=std::count(begin(),begin()+beg,'(')-std::count(begin(),begin()+beg,')');
    if(not open)return "";

    // find start of "(" that encloses term
    size_t left=beg;
    for(;open>0;left--)
        if(std::string(*this)[left]=='(')open--;

    // find end of brackets of that enclose term
    size_t right=left+1;
    for(open=1;open>0;right++){
        if(std::string(*this)[right]=='(')open++;
        if(std::string(*this)[right]==')')open--;
    }
    return substr(left,left-right+1);

}

OperatorDefinition OperatorDefinition::dropTerms(string Hierarchy) const{
    if(*this=="")return *this; // nothing to be dropped
    if(Hierarchy.find("&")!=std::string::npos and not isExpanded())return *this; // do not drop on unresolved hybrid hierarchies

    vector<OperatorDefinition>term=OperatorDefinition(std::string(*this)).terms();
    OperatorDefinition::dropTerms(term,Index::coordinates(Hierarchy));
    string s;
    for(size_t k=0;k<term.size();k++)s+=term[k];
    return s;
}

std::string OperatorDefinition::inputCoordinates(string Hierarchy){
    // remove the LOWER of any duplicate coordinate (i.e. remove the levels in the floor)
    vector<string> iH=tools::splitString(Hierarchy,'.');

    string s;
    for(size_t k=0;k<iH.size();k++){
        std::string thisAx=iH[k];
        if(thisAx[0]=='k')thisAx="spec"+thisAx.substr(1); // treat spec/k pair as FE
        if(std::find(iH.begin(),iH.begin()+k,thisAx)==(iH.begin()+k))s+="."+iH[k];
    }
    return s.substr(1); // remove leading "."

}

std::string OperatorDefinition::permuteFEM(string Hierarchy){
    //    vector<string>term=terms(*this);
    vector<OperatorDefinition>term=OperatorDefinition(std::string(*this)).terms();
    string s;
    for(size_t k=0;k<term.size();k++){
        s+=OperatorDefinition(term[k]).axisToCoor(Hierarchy);
    }
    return s;
}

std::string OperatorDefinition::permute(string Term, string InCoor, string OutCoor){
    if(InCoor==OutCoor)return *this;

    vector<string>inp=tools::splitString(InCoor,'.');
    vector<string>out=tools::splitString(OutCoor,'.');
    if(inp.size()!=out.size())ABORT("numbers of coordinates differ: "+InCoor+" --- "+OutCoor);
    vector<unsigned int>perm;
    for(int k=0;k<inp.size();k++)
        perm.push_back(std::find(inp.begin(),inp.end(),out[k])-inp.begin());
    if(*std::max_element(perm.begin(),perm.end())==inp.size())
        ABORT("coordinates not related by permutation: "+InCoor+" --- "+OutCoor);
    return opPermuted(Term,perm);
}

std::string OperatorDefinition::compose(string Factor, string Term){
    if(Term=="" or Factor=="")return Factor+Term;
    if(Term.find("$")!=string::npos){
        // no explict constraint on factor, assume none
        if(Factor.find("$")==string::npos)
            Term.insert(Term.find(".",Term.find("$")),".*");
        else {
            // move factor constraint to overall
            Term.insert(Term.find(".",Term.find("$")),Factor.substr(Factor.find(".",Factor.find("$"))));
            Factor=Factor.substr(Factor.find("$"));
        }
    }
    else if(Factor.find("$")!=string::npos)
        ABORT(Str("if factor has a constraint, term must as well, found:")+Factor+Term);

    return Factor+Term;
}
