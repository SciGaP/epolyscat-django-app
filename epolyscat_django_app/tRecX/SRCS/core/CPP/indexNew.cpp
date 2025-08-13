// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "indexNew.h"

#include "basisSub.h"
#include "readInput.h"
#include "printOutput.h"
#include "overlapDVR.h"
#include "inverseDvr.h"
#include "inverseFem.h"
#include "axisTree.h"
#include "basisOrbital.h"
#include "basisIntegrableSqrt.h"
#include "basisIntegrableRpow.h"
#include "indexConstraint.h"
#include "basisDvr.h"
#include "basisVector.h"
#include "basisExpIm.h"
#include "indexOverlap.h"
#include "basisAbstract.h"

#include "eigenTools.h"
#include "algebra.h"

using namespace std;

bool IndexNew::doNotUse=false;
static std::map<const IndexNew*,const AxisTree*> _axisTree;
IndexNew::~IndexNew(){
    if(not parent() and _axisTree.count(this))
        _axisTree.erase(this);
}
/// resolve various dependencies
void IndexNew::resolveDependence(BasisSetDef & Def,std::vector<unsigned int> Branch, const std::vector<const Index*> Path){
    BasisIntegrableSqrt::resolvePower(Def,Branch,Path);
    BasisIntegrableRpow::resolvePower(Def,Branch,Path);
}

static std::vector<BasisSetDef> modifyFE(const std::vector<BasisSetDef>& BasDef, std::vector<const Index*> Path){
    std::vector<BasisSetDef> res;
    std::vector<std::string> keywords({"order","only","ommit","keep"});
    for(auto b: BasDef){
        bool accept=true;
        std::vector<std::string> mods=tools::splitString(b.modify(),'&');
        // loop through all modifiers
        for(auto mod: mods){
            mod=tools::cropString(mod);
            for(auto k: keywords)
                if(mod.find(k)!=std::string::npos and mod.find(k)!=0)
                    ABORT("found second modifier \""+k+"\" in \""+mod+"\", separate by & (for logical \"and\")");
            if(mod.find("ommit")==0){
                accept=false;
            }
            else if(mod.find("only ")==0){
                // axis index
                std::vector<std::string> axCond=tools::splitString(tools::cropString(mod.substr(mod.find("only ")+5)),' ');
                if(axCond.size()!=2)ABORT("need format 'only ax1.ax2..axN condition(s)], got "+mod);
                auto pdep=std::find_if(Path.begin(),Path.end(),[&](const Index* Idx){return Idx->axisName()==axCond.front();});
                if(pdep==Path.end())ABORT("cannot find axis "+axCond.front()+" for condition "+mod);
                auto dep=*pdep;
                if(tools::cropString(axCond[1]).find("L<")==0){
                    Algebra lim(tools::cropString(axCond[1]).substr(2));
                    if(not lim.isAlgebraOfConsts())ABORT("not a well-formed constant expression after L< in "+mod);
                    // current branch will be appended next to dep
                    if(dep->basis()->physical(dep->childSize())>=lim.val(0.).real())accept=false;

                }
                else ABORT("cannot interprete condition in "+mod);
            }
            else if(mod.find("order=")==0){
                std::vector<std::string> axCond=tools::splitString(tools::stringInBetween(mod,"[","]",true),':');
                if(axCond.size()==1 and axCond[0]!=""){
                    b.order=tools::string_to_int(mod.substr(mod.find("=")+1));
                }
                else ABORT("cannot interprete condition in "+mod);
            }
            else if (mod.find("keep")!=0)
                DEVABORT("cannot interpret modification string "+mod);

        }
        if(accept)res.push_back(b);
    }

    // (re-)establish first/last indicators (little clumsy)
    // case: multiple elements, intermediate elements removed
    for(size_t k=1;k<res.size();k++){
        if(res[k-1].upBound()!=res[k].lowBound()){
            res[k-1].last=true;
            res[k].first=true;
        }
    }
    // case: multiple elements, first element removed, impose boundary as on old first
    if(res.size() and res.front().lowBound()!=BasDef.front().lowBound())res.front().first=BasDef.front().first;

    // case: multiple elements, last element removed, impose boundary as on old last
    if(res.size() and res.back().upBound()!=BasDef.back().upBound())res.back().last=BasDef.back().last;
    return res;
}

void IndexNew::addNonEmpty(AxisTree* Ax, const IndexConstraint* Constraint,std::vector<unsigned int>  Pos,
                           std::vector<const Index *> Path,std::vector<int> & Removed)
{
    if(Constraint !=0 and not Constraint->includes(Path, Pos)){
        Removed.push_back(Pos.back());
    }
    else {
        Index* child(0);
        try{child=new IndexNew(Ax,Constraint,Pos,Path,false);
        }catch(Index::empty_subtree_exception& ex){
            delete child;
            child=0;
        }
        if(child != 0){
            if(not noDum or Ax)
                childAdd(child);
            else {
                _size=basis()->size();
                delete child;
            }

        }else{
            Removed.push_back(Pos.back());
            delete child;
        }
    }

}

void IndexNew::setBasisRemoved(const std::vector<int> & Removed){
    if(not basis() or Removed.size() == basis()->size()){
        throw Index::empty_subtree_exception();
    }
    if(Removed.size()>0 and Removed.size()!=basis()->size()){
        if(dynamic_cast<const BasisDVR*>(basis()))ABORT("cannot constrain DVR basis (for now)");
        setBasis(basis()->remove(Removed));
    }
}

void IndexNew::buildOverlap() {
    // compute overlap matrices in places where it makes sense
    localOverlapAndInverse(0,0);
    IndexOverlap::set(this);
    setInverseOverlap(Inverse::factory(this));
}


IndexNew::IndexNew(const AxisTree *Ax, const IndexConstraint* Constraint,
                   std::vector<unsigned int> InPos, std::vector<const Index *> InPath, bool Entry)
    :IndexNew()
{
    std::vector<unsigned int> Pos(InPos);
    std::vector<const Index *> Path(InPath);
    if(Entry){
        _axisTree[this]=Ax;
    }
    std::vector<int> removed;
    if(Ax==0){
        // end of axis hierarchy
        setBasis(BasisAbstract::factory("Vector:1"));
        setAxisName("NONE");
        unsetFloor();
    }
    else {
        std::vector<BasisSetDef> basDef=modifyFE(Ax->basDef,Path);

        Pos.push_back(0);
        Path.push_back(this);
        setAxisName(Ax->name);
        setAxisSubset(Ax->subset());

        if(Ax->childSize()>1){
            // hybrid axis
            if     (Ax->bases.size())setBasis(Ax->bases[0]);
            else if(basDef.size()){
                setBasis(BasisAbstract::factory(basDef[0]));
            }//OBSOLESCENT style

            Index* hybIdx=this;
            if(Ax->name.find("&")==std::string::npos){
                // if not explicit hybrid axis, need extra Index level for hybrid
                setBasis(BasisAbstract::factory("Hybrid: 1"));
                childAdd(new Index());
                hybIdx=child(0);
                Path.push_back(hybIdx);
                Pos.push_back(0);
            }
            string hybsub=Ax->child(0)->subset();
            string hybnam=Ax->child(0)->name;
            for(size_t k=1;k<Ax->childSize();k++){
                hybnam+="&"+Ax->child(k)->name;
                hybsub+="&"+Ax->child(k)->subset();
            }
            hybIdx->setAxisName(hybnam);
            hybIdx->setAxisSubset(hybsub);
            hybIdx->setBasis(BasisAbstract::factory(Sstr+"Hybrid:"+Ax->childSize()));
            for(size_t k=0;k<Ax->childSize();k++){
                Pos.back()=k;
                hybIdx->childAdd(new IndexNew(Ax->child(k),Constraint,Pos,Path,false));
            }
            if(hybIdx!=this){
                Pos.pop_back();
                Path.pop_back();
            }
        }

        else if(basDef.size()==1){
            // standard, non FE axis
            if(Ax->bases.size()==1)
                setBasis(Ax->bases[0]);
            else {
                BasisSetDef curDef=basDef[0].resolveDependence(Pos,Path); // OBSOLESCENT - do not add functionality
                resolveDependence(curDef,Pos,Path); // new form - move over from BasisSetDef
                setBasis(BasisAbstract::factory(curDef));
            }
            if(not basis() or basis()->size()==0)DEVABORT(Sstr+"empty basis"+basis());
            if(basis()->size()==0)ABORT("zero-size basis on "+strNode());

            // descend to next level axes
            for(unsigned int k=0;k<basis()->size();k++){
                Pos.back()=k;
                addNonEmpty(Ax->descend(),Constraint,Pos,Path,removed);
            }
        }
        else if(basDef.size()>1){
            // finite element axis
            std::vector<int> subs;

            for(size_t k=0; k<Ax->basDef.size();k++){
                if(subs.size()<basDef.size() and
                        Ax->basDef[k].lowBound()==basDef[subs.size()].lowBound() and
                        Ax->basDef[k].upBound()==basDef[subs.size()].upBound())subs.push_back(k);
            }
            if(subs.size()!=basDef.size())DEVABORT("wrong subset");

            // we may add a BasisFE... to replace Dummy
            auto bas=BasisAbstract::factory("Vector: +"+tools::str(Ax->basDef.size()));
            setBasis(bas->size()==subs.size()?bas:BasisAbstract::factory(BasisSub::strDefinition(bas,subs)));
            if(not basis() or basis()->size()==0)DEVABORT(Sstr+"empty basis"+basis());
            AxisTree * axt=const_cast<AxisTree*>(Ax);
            for(unsigned int k=0;k<basis()->size();k++){
                // append single element floor axis
                Pos.back()=k;
                BasisSetDef bDef=basDef[k];
                resolveDependence(bDef,Pos,Path); // new form - move over from BasisSetDef

                //not very pretty...
                axt->firstLeaf()->childAdd(new AxisTree(Axis(Ax->name,Ax->comsca,bDef)));

                childAdd(new IndexNew(axt->child(0),Constraint,Pos,Path,false));
                axt->descend(axt->height()-1)->childPop();
            }
            // remove possible empty subtrees
            for(size_t k=0;k<childSize();k++){
                if(not child(k)->basis())removed.push_back(k);
            }
        }
        setBasisRemoved(removed);
        Pos.pop_back();
        Path.pop_back();
    }


    if(Path.size()==0){
        sizeCompute();
        std::vector<std::string> feAx;
        setFloorAuto(feAx);

        // register top node of each axisSubset
        for (const Index* idx=this;idx!=0;idx=idx->nodeNext()){
            if(not BasisOrbital::getReference(idx->axisSubset())){
                BasisOrbital::addIndex(idx->axisSubset(),idx);
            }
        }
        buildOverlap();
    }
}

const AxisTree* IndexNew::axes() const{
    if(not _axisTree.count(this))return 0;
    return _axisTree[this];
}

void IndexNew::read(ReadInput& Inp){

    std::vector<std::string> axisCategories=Inp.categoryWithSpecification("Axis");
    axisCategories.push_back("Axis");

    for(auto ax: axisCategories){
        // loopy through Axis and possibly Axis[...]
        if(tools::cropString(ax)=="")continue;

        // seek and set up alternate versions of main Axis
        std::set<std::string> alternates=Inp.namesWithSpecification(ax,"alternateAxis");
        alternates.insert("alternate");
        for(auto a: alternates){
            if (tools::cropString(a)!="") {

                // advance to first line not to be omitted and with non-empty alternative
                int line(0);
                std::string alt("");
                while(not ReadInput::main.endCategory("Axis",++line)){
                    alt=Axis::readAlternate(ReadInput::main,line,ax,"keep",a);
                    if(alt!="ommit")break;
                }
                if(alt!=""){
                    // alternate was actually specified - build and register index
                    AxisTree axt(Inp,line,"",ax,a);
                    BasisOrbital::addIndex(a,new IndexNew(&axt));
                    PrintOutput::message("added \""+a+"\" to list of reference indices");
                }
            }
        }
    }

    //    for(auto n: Inp.categoryWithSpecification("Axis")){
    //        if(tools::cropString(n)!=""){
    //            int line(1);
    //            AxisTree ax(Inp,line,"",n);
    //            BasisOrbital::addIndex(tools::stringInBetween(n,"[","]"),new IndexNew(&ax));
    //            PrintOutput::message("added \""+tools::stringInBetween(n,"[","]")+"\" to list of reference indices");
    //        }
    //    }
}


