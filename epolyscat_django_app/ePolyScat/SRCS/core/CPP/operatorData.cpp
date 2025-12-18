// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorData.h"
#include "readInput.h"
#include "stringTools.h"
#include "abort.h"
#include "useMatrix.h"
//#include "basisMat.h"
#include "basisDvr.h"
#include "parameters.h"
#include "operatorSingle.h"
#ifdef _USE_HACC_
#include "discretizationHaCC.h"
#endif

using namespace tools;
using namespace std;

// Extensions for left and right substitution names
const static string lExt = "left";
const static string rExt = "right";
map<string,map<string,OperatorDefinition> > OperatorData::standardOperators;


void OperatorData::factorize(string Def, Discretization *IDisc, Discretization *JDisc, Index *IIndex, Index *JIndex,
                             std::vector<string> &Factors, std::vector<string> &Remainders){
    // split Def into a first factor(s) and the remainder(s)
    // rules:
    // on any level, an operator is a sum of terms
    // the number of tensor factors must match the number of axes
    // the leading factor must be single-coordinate
    // if leading factor is enclosed in round brackets, split into several terms
    // continuity levels are ignored, the corresponding floor levels must be the last factors in the definition
    // for now, only local operators (multiplication, differentiation) are allowed
    // if no axes are defined, return empty string (= no tensor factor)

    // leave definition unchanged and return factor <Id>
    // if:
    // continuity level
    // left and right levels differ
    // left and right discretizations differ
    // no axes are defined (i.e. non-standard basis)

    if(Def=="mapFromParent" or Def=="mapToParent"){
        Factors.assign(1,Def);
        Remainders.assign(1,Def);
        return;
    }

    if(IIndex->depth()!=JIndex->depth() or
            tools::anyElement(JDisc->continuityLevel,tools::equal,int(IIndex->depth())) or
            IDisc->axis.size()==0){
        Factors.assign(1,"<Id>");
        Remainders.assign(1,Def);
        return;
    }

    // clean up possible white spaces
    Def=tools::cropString(Def);

    // split everything outside <> or () brackets
    vector<string> term0=OperatorData::terms(Def);

    Factors.clear();
    Remainders.clear();
    for (unsigned int n=0;n<term0.size();n++){
        if(Def.find('(')<Def.find('<')){
            // factor is a sum of terms, e.g. (<dJd>+<qJq>)<J>
            if(Def.find('(')!=0)ABORT("no factors outside round brackets allowed: "+Def);
            if(Def.find('(',1)<Def.find(')',1))ABORT("no nested round brackets allowed: "+Def);
            vector<string>facs=OperatorData::terms(term0[n].substr(1,term0[n].find(')')-1));
            string rem=term0[n].substr(term0[n].find(')')+1);
            for (unsigned int k=0;k<facs.size();k++){
                Factors.push_back(facs[k]);
                Remainders.push_back(rem);
            }
        } else {

            Factors.push_back(term0[n].substr(0,term0[n].find('>')+1));
            Remainders.push_back(term0[n].substr(term0[n].find('>')+1));
        }
    }
}

bool OperatorData::isStandard(const string &Def){
    if(Def.find(">")==string::npos) return false;
    vector<string>terms,signs;
    splitString(Def,"+-",terms,signs,"<|",">>");
    return terms.size()!=0 and terms[0].find_first_of("<|")!=string::npos;
}

bool OperatorData::isMultip(const string &Def){
    if(Def.find(">")==string::npos) return false; // not standard form
    if(tools::subStringCount(Def,">")!=2) return false; // for now, only 2d
    if(Def.find("d_")!=string::npos or Def.find("_d")!=string::npos) return false; // not for derivatives (although...)
    return true;
}

bool OperatorData::isZero(const string &Def, const Index *IIndex, const Index *JIndex){
    if(Def.find("ZERO")!=string::npos)return true;
    return false;
}

vector<string> OperatorData::singleTerms(string Def,const std::string & IHierarchy, const std::string & JHierarchy){
    // unbracket and split
    vector<string> term=terms(OperatorDefinition::unBracket(Def).str());
    string unbra;
    for(int n=0;n<term.size();n++){
        // expand if first term is predefined operator
        if(term[n].find("<<")==term[n].find("<")){
            unbra+=expandStandard(term[n],IHierarchy,JHierarchy);
        }
        else {
            unbra+=term[n];
        }
    }
    // expansion may have created new brackets
    return terms(OperatorDefinition::unBracket(unbra).str());
}

// split at "+-" outside () or <>
vector<string> OperatorData::terms(string Def){
    vector<string> term,sign,ret0,ret;
    Def=cropString(Def);

    if(Def.substr(0,1)=="'"){
        if(Def.find("'",1)!=Def.length()-1)ABORT("unmatched first ' in "+Def);
        Def=Def.substr(1);
        Def=Def.substr(0,Def.length()-1);
    }

    // strip overall round brackets
    if(Def[0]=='(' and Def[Def.length()-1]==')')Def=Def.substr(1,Def.length()-2);

    splitString(Def,"+-",term,sign,"<|(",">>)");
    for(unsigned int n=0;n<term.size();n++)
        if(term[n].length()!=0 and term[n]!=" ")ret0.push_back(cropString(sign[n]+term[n]));

    // un-bracket grouped terms
    for(unsigned int k=0;k<ret0.size();k++){
        if(ret0[k].rfind(")")==ret0[k].length()-1){
            string par=cropString(ret0[k].substr(0,ret0[k].find("(")));
            splitString(stringInBetween(ret0[k],"(",")"),"+-",term,sign,"<|(",">>)");

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
            ABORT("for now, cannot have nested (...) in operator: "+Def);

    return ret;
}

/// \param Term contains depencencis in the form termName{ax1,ax2} etc.
/// \param Dep  with elements ax1,ax2 etc.
void OperatorData::dependence(const string &Term, std::vector<string> &Dep){
    Dep.clear();
    if(Term.find('{')==string::npos)return; // does not have dependencies
    if(subStringCount(Term,"{")>subStringCount(Term,"{}")+1)
        ABORT("cannot have multiple dependencies in single term: "+Term);
    string Ts=Term.substr(Term.rfind("<"));
    Dep=tools::splitString(tools::stringInBetween(Ts,"{","}"),',');
    return;
}

string OperatorData::parameter(const string & Term){

    // make sure basic parameters are defined
    Parameters::defaults();

    // check that SINGLE term
    if(OperatorData::terms(Term).size()>1 and Term.find_first_of("<|")<Term.find("("))return "1";
    size_t start=min(Term.find_first_of("<|("),Term.find("[["));
    if(start==string::npos)return "1";

    // point to overall factor
    string par=cropString(Term.substr(0,start));
    par=par.substr(0,max(max(par.rfind("."),par.rfind("<")),par.rfind("("))); // flexible syntax: allow '.'
    par=par.substr(0,par.find("[[")); // [[...]] indicates "special" operator
    return par;
}

string OperatorData::sign(const string & Def){
    vector<string> elem,seps;
    tools::splitString(Def,"+-",elem,seps);
    return seps[0];
}

string OperatorData::first(const string & Def,bool unsign){
    if(not isStandard(Def))return Def;

    string def=Def;
    if(unsign){
        vector<string>elem,seps;
        tools::splitString(def,"+-",elem,seps,"(<",")>");
        def=elem[0];
    }
    return def.substr(0,def.find(">")+1);
}

string OperatorData::remainder(const string & Def){
    if(not isStandard(Def))return Def;
    return Def.substr(Def.find(">")+1);
}

void OperatorData::substitutionIndices(const string &Def, substitutionSpec& spec, const string& leftName, const string& rightName,
                                       std::vector<size_t> &leftSub, std::vector<size_t> &rightSub)
{
    // clear left-over entries
    leftSub.clear(); rightSub.clear();

    // find leftName
    string temp(leftName);
    temp.insert(0, lExt); // temp stores left[leftName]
    size_t strPos = Def.find(temp);
    while (strPos!=string::npos) {
        leftSub.push_back(strPos);
        leftSub.push_back(temp.size());
        strPos = Def.find(temp, strPos+1);
    }

    // find rightName
    temp = rightName;
    temp.insert(0, rExt); // temp stores right[rightName]
    strPos = Def.find(temp);
    while (strPos!=string::npos) {
        rightSub.push_back(strPos);
        rightSub.push_back(temp.size());
        strPos = Def.find(temp, strPos+1);
    }

    // set substitutionSpec
    if (rightName==string("surfR") or rightName==string("surfX") or rightName==string("surfY") or rightName==string("surfZ")
            or rightName==string("surfRn") or rightName==string("surfRn1") or rightName==string("surfRn2")) {
        spec=surfaceOnFloor;} // information on floor
    else { spec=index; }

    return;
}

void OperatorData::setStandard(){OperatorDefinition::setup();}

string OperatorData::expandStandard(const string Def, const string  & IHierarchy, const string & JHierarchy){

    if(Def.find("<<")==string::npos)return Def; // nothing to be expanded

    string def=OperatorDefinition::unBracket(Def).str();
    if(terms(def).size()!=1)ABORT("only for expanding single term, is: "+Def);
    string coor=Index::coordinates(IHierarchy);
    if(coor!=Index::coordinates(JHierarchy)){
        size_t beg;
        while((beg=def.find("<<"))!=string::npos){
            size_t end=def.find(">>");
            if(end==string::npos)ABORT("unmatched <<...>>: "+def);
            def.replace(beg,2,"[[");
            def.replace(end,2,"]]");
        }
    }
    else {
        def=expandStandard(def,IHierarchy);
    }
    return def;
}

string OperatorData::expandStandard(const string Def, const Discretization *Disc){return expandStandard(Def,Disc->idx()->hierarchy());}
string OperatorData::expandStandard(const string Def, std::string Hierarchy){

    // nothing to be expanded
    if(Def.find("<<")==string::npos)return Def;

    // create the table of standard operators
    setStandard();

    vector<string> term,sep;
    tools::splitString(OperatorDefinition::unBracket(Def).str(),"+-",term,sep,"<(",">)");

    vector<string> axis=tools::splitString(Hierarchy,'.');

    //HACK temporary, until more  general handling
    int kax=0;
    for(;kax<int(axis.size());kax++){
        int l=0;
        for(;l<kax;l++){
            if(axis[l]==axis[kax]){
                break; // repeated appearance of axis indicates floor
            }
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
    if(ReadInput::main.flag("DEBUGshowOp","show the expanded operator definitions"))PrintOutput::lineItem(Def+" expands to: ",def);
    return def;
}

void OperatorData::constrain(UseMatrix &Mult, string Constr, unsigned int Level){
    size_t cPos=Constr.find('$');
    if(cPos==string::npos)return;

    vector<string>lConst=tools::splitString(Constr.substr(cPos),'.');
    if(lConst.size()<Level+2)ABORT("number of constraints does not match level="+tools::str(Level)+", def="+Constr);

    if(Constr.substr(cPos).find("$BAND")==0){
        if(lConst[Level+1]=="*")return;
        if(string("0123456789").find(lConst[Level+1])==string::npos)ABORT("illegal band width specification at level="+tools::str(Level)+" in "+Constr);
        int c=string_to_int(lConst[Level+1]);
        for(int j=0;j<Mult.cols();j++)
            for(int i=0;i<Mult.rows();i++)
                if(abs(i-j)>c)Mult(i,j)=0.;
        unsigned int truesub,truesuper,nonz;
        char kind;
        Mult.diagnose(0.,nonz,truesub,truesuper,kind);
        if(int(truesub+truesuper)>2*c)cout<<"constrain failed: "<<lConst[Level+1]<<" c="<<c<<" multiplier band: "<<truesub+truesuper<<endl;
    }
    else
        ABORT("undefined constraint in definition "+Constr);
}

string OperatorData::extractBlock(string Definition, unsigned int I, unsigned int J){
    string extract;
    vector<string> trms(terms(Definition));
    for(int k=0;k<trms.size();k++)
    {

        // adjust constraints (if any)
        if(trms[k].find("$")!=string::npos){
            // first constraint parameter and remove
            size_t pos0=trms[k].find(".",trms[k].find("$"));
            size_t pos1=min(trms[k].length(),trms[k].find(".",pos0+1));
            trms[k].erase(pos0,pos1-pos0);
        }

        // loop through terms
        string blk=first(trms[k]); // first factor of term
        size_t posL=blk.find("<");
        size_t posR=blk.find(">");
        size_t posM=blk.find(",");
        if(posM==string::npos and blk.find("<1>")==string::npos and  blk.find("<allOnes>")==string::npos)
            ABORT("term["+tools::str(k)+"] = "+trms[k]+" of "+Definition+" is not block-index: "
                  +blk+" must of the form <i,j>, <1>, or <allOnes>");

        // if matches <I,J> (or <1> or <allOnes>), remove first factor and append to block definition
        if(blk.find("<allOnes>")==posL or (I==J and blk.find("<1>")==posL))
            extract+=trms[k].substr(0,posL)+trms[k].substr(posR+1);
        else {
            auto si = blk.substr(posL+1,posM-posL-1);
            auto sj = blk.substr(posM+1,posR-posM-1);
            auto ci = tools::string_to_int(si);
            auto cj = tools::string_to_int(sj);
            if(posM != string::npos && ci ==I and cj==J){
                extract+=trms[k].substr(0,posL)+trms[k].substr(posR+1);
            }
        }

    }
    return extract;
}

void OperatorData::operatorData(string Def, const Index * const IFloor, const Index * const JFloor, vector<UseMatrix> & Mats)
{

    if(Def.find("NONORTH")!=string::npos){
#ifdef _USE_HACC_
        DEVABORT("DiscretizationHaCC::operatorData no longer available")
//        //HACK for compatibility with old DiscretizationHaCC
//        DiscretizationHaCC::current->operatorData(Def,IFloor,JFloor,Mats);
        return;
#else
        DEVABORT("for haCC, compile with -D_USE_HACC_");
#endif
    }

    const Index * iIndex=IFloor->root();

    Mats.clear();
    vector<complex<double>*> dummyPointer;
    vector<string> dep;
    OperatorData::dependence(Def,dep);

    if(dep.size()<2){
        DEVABORT("disabled");
//        OperatorData::get(Def,IFloor->basProd(),JFloor->basProd(),Mats,dummyPointer);
    }
    else if(dep.size()>2) {
        ABORT("at most 2-dimensional floors implemented");
    }
    else if (dep.size()==2){

        Index * i0=iIndex->axisIndex(dep[0]),*i1=iIndex->axisIndex(dep[1]);
        if(i0!=0 and i0->depthOfDuplicate()!=Index::npos and i1!=0 and i1->depthOfDuplicate()!=Index::npos){
            // 2-dimensional floor
            if(not BasisDVR::femDVR)ABORT("non-tensor product on floor level only for femDVR");
            DEVABORT("disabled");
//            OperatorData::get(Def,IFloor->basProd(),JFloor->basProd(),Mats,dummyPointer);
        }
        else {
            ABORT("this path is disabled - reactivate if needed");
            //            vector<const BasisSet*>iBas,jBas;
            //            vector<unsigned int> iSub,jSub;
            //            // get basis functions and index related to dependence
            //            for(unsigned int k=0;k<dep.size();k++){
            //                unsigned int lev=tools::locateElement(hierarchy,dep[k]);
            //                // get basis set
            //                vector<const Index *>iPath(IFloor->path()),jPath(JFloor->path());
            ////                unsigned int floorLev=tools::locateElement(continuityLevel,int(lev));
            //                unsigned int floorLev=IFloor->i;
            //                if(floorLev<continuityLevel.size()){
            //                    iBas.push_back(IFloor->basProd()[floorLev]);
            //                    jBas.push_back(JFloor->basProd()[floorLev]);
            //                } else {
            //                    iBas.push_back(iPath[lev]->basisSet());
            //                    jBas.push_back(jPath[lev]->basisSet());
            //                    iSub.push_back(iPath[lev]->nSub(iPath[lev+1]));
            //                    jSub.push_back(jPath[lev]->nSub(jPath[lev+1]));
            //                }
            //            }
            //            OperatorData::get(def,iBas,jBas,Mats,dummyPointer,iSub,jSub);
        }
    }
}


/// for all finite element axes, the elements form the floors, while element indices remain on current level <br>
/// the function modifies the definition by moving the factors from the finite element axis level to floor<br>
/// on the current level, an identity matrix is inserted<br>
/// Example: <J><dJd><Jd> ==> <J><Id><Jd><dJd> if the second axis is finite elements
void OperatorData::modifyDefinition(string &Def, const Index* IIndex, const Index* JIndex){

    if(IIndex->heightAboveBottom()==0 and IIndex->depth()==0) return;  // 1D discretization no need to modify

    // WARNING: Code doesn't work with a mixed floor - (with and without continuity)

    // make Def to match with index hierarchy (insert operation on continuity level)
    if(not isStandard(Def))return;
    if(IIndex->firstFloor()->depth()!=JIndex->firstFloor()->depth())
        ABORT("discretization depth do not match: "+IIndex->hierarchy()+" vs. "+JIndex->hierarchy()+"\n cannot use standared operator definition: "+Def);

    string inDef=Def;
    string constraint;
    if(Def.find('$')!=string::npos){
        constraint=Def.substr(Def.find('$'));
        Def=Def.substr(0,Def.find('$'));
    }

    // split into terms
    vector<string> terms=OperatorData::terms(Def);
    //for (int i=0; i!=terms.size(); ++i) { cout << terms.at(i) << endl; }
    Def="";

    for (unsigned int n=0;n<terms.size();n++){
        Index *iI=const_cast<Index*>(IIndex),*jI=const_cast<Index*>(JIndex);
        string mDef=OperatorData::sign(terms[n]),fDef="";
        while((not iI->hasFloor()) and (not jI->hasFloor())){
            //            if(iDisc->isContinuous(iI->depth())){
            if(iI->depthOfDuplicate()!=Index::npos){
                //                if(not jDisc->isContinuous(iI->depth()))ABORT("continuity levels do not match");
                if(iI->depthOfDuplicate()!=jI->depthOfDuplicate())ABORT("continuity levels do not match");
                // note: for now we assume local operator on continuity level
                if(inDef.find("NONORTH")==string::npos)
                    mDef.append("<Id>");                              // continuity level: unity
                else
                    mDef.append("<allOnes>");                         // continuity level: nonorthogonality in haCC mixes finite elements
                fDef.append(OperatorData::first(terms[n],true)); // attach first term stripped of its sign
            }
            else if((iI->axisName().substr(0,4)=="surf" and jI->axisName().substr(0,4)=="surf") or
                    (iI->axisName().substr(0,4)=="spec" and jI->axisName().substr(0,4)=="spec")){
                mDef.append("<Id>");                              // continuity level: unity
                fDef.append(OperatorData::first(terms[n],true)); // attach first term stripped of its sign
            }
            else {
                mDef.append(OperatorData::first(terms[n],true));
            }
            terms[n]=OperatorData::remainder(terms[n]);
            iI=iI->child(0);
            jI=jI->child(0);
        }

        //HACK
        if (fDef.empty()) { // no continuity level, i.e. the floor was given in the definition
            //            fDef.append(OperatorData::first((terms[n])));
            fDef=tools::cropString(terms[n]);
            if(fDef[0]=='+' or fDef[0]=='-')fDef=fDef.substr(1);
        }

        Def.append(mDef.append(fDef));
    }

    Def+=constraint;
}
