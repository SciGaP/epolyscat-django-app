// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "discretization.h"

//resolve forward declarations
#include "wavefunction.h"
#include "index.h"
#include "multiIndex.h"
#include "axis.h"
#include "printOutput.h"
#include "parallel.h"
#include "operatorData.h"
#include "discretizationHybrid.h"
//#include "basisMat.h"
#include "indexConstraint.h"
#include "visualizeAngularDistribution.h"

#ifdef _USE_HACC_
#include "discretizationHaCC.h"
#endif

#include <memory>
#include "memInfo.h"
#include "operatorFloor.h"
#include "axisTree.h"
#include "indexNew.h"
#include "overlapDVR.h"
#include "basicDisc.h"
#include "basisDvr.h"

using namespace std;
using namespace tools;

Discretization::~Discretization(){
}

Discretization * Discretization::factory(ReadInput &Inp){

    if(DiscretizationHybrid::isHybrid(Inp))
        return new DiscretizationHybrid(Inp);
    vector<Axis> ax;
    Axis::fromFile(Inp,ax);
//    for(size_t k=0;k<ax.size();k++)
//        if(ax[k].name=="Ion"){
//#ifdef _USE_HACC_
//            return new DiscretizationHaCC(Inp);
//#else
//            DEVABORT("compile with -D_USE_HACC_");
//#endif
//        }
    IndexConstraint::main.read(Inp);

    return new BasicDisc(Inp,"", &IndexConstraint::main);
}



string Discretization::str(unsigned int Brief) const{

    string s,b;
    s="Discretzation: "+name+"\n";
    if(parent!=0)s+="Parent discretization: "+parent->name+"\n";
    s+="Hierarchy: "+hierarchy[0];
    b=hierarchy[0];
    for(unsigned int n=1;n<hierarchy.size();n++){
        if(Brief==0)s+="->"+hierarchy[n];
        if(n>0)b+=".";
        b+=hierarchy[n];
    }
    s+="\n";
    s+="Continuity levels: ";
    for(unsigned int d=0;d<continuityLevel.size();d++)s+=" "+tools::str(hierarchy[continuityLevel[d]]);
    s+="\n";
    b+="=";
    for(unsigned int n=0;n<axis.size();n++){
        s+=axis[n].str(Brief)+"\n";
        if(n>0)b+=".";
        b+=axis[n].str(Brief);
    }

    if(idx()!=0)s+=idx()->str();

    if(Brief>0)return b;
    return s;

}

void Discretization::print(string File, string Title) const {

    if(Title==""){
        if(BasisDVR::femDVR)Title="DVR";
        else                Title="FEM";
        Title+=" --- DISCRETIZATION";
    }
    PrintOutput::title(Title);
    PrintOutput::paragraph();

    PrintOutput::lineItem("Name",name);
    if(parent!=0)PrintOutput::lineItem("Parent",parent->name);

    PrintOutput::end();

    PrintOutput::paragraph();
    if(axis.size()>0)Axis::print(axis);
    else {
        PrintOutput::lineItem("Coordinates",idx()->coordinates());
        PrintOutput::paragraph();
    }


    if(constraint!=0){
        PrintOutput::paragraph();
        if(constraint->constraints().size()==1)
            PrintOutput::lineItem("Constraint",constraint->constraints()[0]->str());
        else if (constraint->constraints().size()>1){
            PrintOutput::newRow();
            PrintOutput::rowItem("Constraints");
            for(int k=0;k<constraint->constraints().size();k++){
                PrintOutput::newRow();
                PrintOutput::rowItem(constraint->constraints()[k]->str());
                PrintOutput::rowItem(constraint->constraints()[k]->str());
            }
        }
    }
}

Coefficients& Discretization::inverseOverlap(Coefficients &coeffs,Index* index, bool &CoeffsPerformInvOv) {
    CoeffsPerformInvOv = true;
    return coeffs;
}

void Discretization::setAxis(ReadInput &In, string Subset){
    // get the axes
    Axis::fromFile(In,axis,Subset); // read all axes

    if(axis.size()==0){
        string mess="no axis found";
        if(Subset!="")mess+=" for subset \""+Subset+"\"";
        ABORT(mess);
    }


    // set the hierarchy names
    for(unsigned int n=0;n<axis.size();n++){
        hierarchy.push_back(axis[n].coor.name());
        if((axis[n].basDef.size()>1 and axis[n].basDef[0].name()=="useIndex"))DEVABORT("\"useIndex\" is obsolet, use \"vector\" instead")
        if((axis[n].basDef.size()>1 and axis[n].basDef[0].name()!="vector") or
                (axis[n].basDef.size()==1 and axis[n].basDef[0].name()=="grid"))
            continuityLevel.push_back(n); // multiple elements require continuity
        if(n>0)name+=".";
        name+=axis[n].coor.name();
    }
    if(continuityLevel.size()==0){
        if(axis.back().name.find("Neutral")==string::npos and
                axis.back().name.find("Ion")==string::npos
                )
        {
            PrintOutput::warning("did not find continuity level - assume last level as single element axis");
        }

    }
}

/// generate name from axis names, set up indices, basis, inverse etc.
void Discretization::construct(){
    // all split coordinates go on floor level
    hierarchy.push_back("f");
    for(unsigned int n=0;n<axis.size();n++){
        if(axis[n].basDef.size()>1)hierarchy.back()+=hierarchy[n];
    }
    if(_axisTree==0 and axis.size()>0)_axisTree.reset(new AxisTree(axis));

    Index::build=true;
    idx()=new IndexNew(_axisTree.get(),constraint);

    Inverse::factory(idx());
    Parallel::setSort(this);
    idx()->testInverseOverlap();
    name=idx()->coordinates();
}

string Discretization::coordinates() const {
    return idx()->coordinates();
}
