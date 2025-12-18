// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "indexFromGrid.h"
#include "basisGrid.h"

using namespace std;

IndexFromGrid::IndexFromGrid(const Index *IdxGrid, std::vector<string> Deflate, vector<const Index*> Path)
{

    setAxisName(IdxGrid->axisName());

    Deflate.resize(Path.size(),"1"); // supplement deflate factor where missing
    setBasis(basisFromGrid(IdxGrid,Deflate[IdxGrid->depth()],Path));

    Path.push_back(this);
    for(size_t k=0;k<basis()->size();k++)
        childAdd(new IndexFromGrid(IdxGrid->child(k),Deflate,Path));
}

const BasisAbstract * IndexFromGrid::basisFromGrid(const Index* GridIdx, string ContractFactor, std::vector<const Index*> Path) {
    if(not GridIdx->basis()->isGrid())return GridIdx->basis();
//    const BasisGrid *g=BasisGrid::factory(GridIdx->basis());

    Coordinate coor=Coordinate::fromString(GridIdx->axisName());
    if(coor.defaultFunction()=="useIndex")DEVABORT("\"useIndex\" is obsolete, replace by \"vector\"");
    if(coor.defaultFunction()=="vector")return GridIdx->basis();
    if(not GridIdx->subEquivalent())return GridIdx->basis();

    Str db(coor.defaultFunction(),"");
    vector<double>par;
    if(tools::findFirstOutsideBrackets(db,"{","[","]")){
        string depAx=tools::stringInBetween(db,"{","}"); // dependent axis
        if(GridIdx->axisName().find_first_of("0123456789")!=string::npos)
            depAx+=GridIdx->axisName().substr(GridIdx->axisName().find_first_of("0123456789")); // append numbering (if any)
        for(const Index* depIdx: Path){
            if(depIdx->axisName()==depAx){
                db=Str(db.substr(0,db.find("{")),"")+"["+depIdx->basis()->physical(depIdx->childSize())+"]";
                par.push_back(depIdx->basis()->physical(depIdx->childSize()));
            }
        }
    }
    int order=GridIdx->basis()->size()/tools::string_to_int(ContractFactor);
    BasisSetDef def(order,coor.qmin,coor.qmax-coor.qmin,coor.defaultFunction(),true,true,true,
                    coor,{0,order-1},ComplexScaling(),false,par);

    return BasisAbstract::factory(def);
}
