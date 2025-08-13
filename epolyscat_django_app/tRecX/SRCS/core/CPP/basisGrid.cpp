// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisGrid.h"

#include "basisAbstract.h"
#include "basisDvr.h"
#include "basisSub.h"
#include "interpolate.h"

std::map<std::vector<double>,const BasisGrid*>BasisGrid::_allBasisGrid;

BasisGrid::BasisGrid(const BasisAbstract *Grid)
    :BasisAbstract("grid")
{
    if(not BasisSub::superBas(Grid)->isGrid() and not Grid->isDVR())DEVABORT("input BasisAbstract does not seem to be grid-type: "+Grid->str());

    const BasisGrid * g;
    const BasisDVR * d;
    const BasisSub * s;
    if((g=dynamic_cast<const BasisGrid*>(Grid))!=0)
        _mesh=g->mesh();
    else if((d=dynamic_cast<const BasisDVR*>(Grid))!=0){
        std::vector<double>dum;
        d->dvrRule(_mesh,dum);
    }
    else if((s=dynamic_cast<const BasisSub*>(Grid)) and BasisSub::superBas(s)->isGrid()){
        BasisGrid super(BasisSub::superBas(s));
        for(size_t k=0;k<s->size();k++)_mesh.push_back(super.mesh()[s->subset()[k]]);
    }
    else {
        for (size_t k=0;k<Grid->grid()->size();k++){
            _mesh.push_back(Grid->grid()->mesh()[k]);
        }
    }
}

bool BasisGrid::operator==(const BasisAbstract& Other) const {
    if(not Other.isGrid())return false;

    const BasisGrid * o=BasisGrid::factory(&Other);
    if(o==0)return false;
    if(o->size()!=size())return false;
    double eps;
    eps=std::max(abs(mesh().back()-mesh()[0])*1.e-1,1.e-12);
    for(size_t k=0;k<mesh().size();k++)
        if(abs(mesh()[k]-o->mesh()[k])>eps)return false;
    return true;
}

const BasisGrid* BasisGrid::factory(std::vector<double> Mesh){
    const BasisGrid* g=_allBasisGrid[Mesh];
    if(g==0)g=_allBasisGrid[Mesh]=new BasisGrid(Mesh);
    return g;
}

const BasisGrid* BasisGrid::factory(const BasisAbstract *Grid){
    if(not BasisSub::superBas(Grid)->isGrid())return 0;
    BasisGrid g(Grid);
    return factory(g.mesh());
}

const BasisGrid* BasisGrid::factory(const std::string Definition){
   if(Definition.substr(0,Definition.find(":"))!="Grid")
           ABORT("BasisGrid definition must start as \"Grid:...\", is: "+Definition.substr(0,12)+"....");
   std::vector<std::string> strMesh=tools::splitString(Definition.substr(Definition.find(":")+1),',');
   std::vector<double> mesh;
   for(auto s: strMesh)mesh.push_back(tools::string_to_double(s));
   return factory(mesh);
}
std::string BasisGrid::str(int Level) const{
    std::string s("Grid ");
    s+="["+tools::str(_mesh[0],3)+","+tools::str(_mesh.back(),3)+"] "+tools::str(size())+"["+tools::str(size())+"]";
    return s;
}

Eigen::MatrixXcd BasisGrid::mapInterpolate(const BasisAbstract *Bas, int Order) const{
    const BasisGrid* gBas=dynamic_cast<const BasisGrid*>(Bas);
    if(not gBas)ABORT("cannot interpolate for non-grid basis: "+Bas->str());
    for(double g: gBas->mesh())
        if(g<gBas->mesh()[0] or gBas->mesh().back()<g)
            DEVABORT("grid point not within original grid, cannot interpolate");
    std::vector<std::complex<double>> val;
    Eigen::MatrixXcd res(gBas->size(),size());
    std::complex<double> dumRes;
    for(size_t k=0;k<mesh().size();k++){
        val.assign(mesh().size(),std::complex<double>(0.,0.));
        val[k]=1.;
        InterpolateNewton<double,std::complex<double>> intpol(mesh(),val,Order);
        for(size_t l=0;l<gBas->size();l++)
            res(l,k)=intpol.val(gBas->mesh()[l],dumRes);
    }
    return res;
}


std::string BasisGrid::strDefinition() const {
    std::string s("Grid:");
    for(size_t k=0;k<_mesh.size();k++)s+=tools::str(_mesh[k])+",";
    return s.substr(0,s.size()-1); // trim the last ','
}

