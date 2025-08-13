// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisGridQuad.h"

std::map<std::vector<double>,const BasisGridQuad*>BasisGridQuad::_allBasis;

BasisGridQuad::BasisGridQuad(const BasisAbstract *Grid)
    :BasisGrid(Grid)
{

    const BasisGridQuad * g=dynamic_cast<const BasisGridQuad*>(Grid);
    if(g!=0){
        _weights=g->weights();
    }
    else {
        if(Grid->gridQuad()->weights().size()!=_mesh.size())
            DEVABORT(Str("number of ")+Grid->gridQuad()->weights().size()+"weights does not match grid size of"+_mesh.size());
        for(size_t k=0;k<Grid->size();k++)_weights.push_back(Grid->gridQuad()->weights()[k]);
    }
}

const BasisGridQuad* BasisGridQuad::append(const BasisAbstract *Grid) const{
    const BasisGridQuad * g=dynamic_cast<const BasisGridQuad*>(Grid);
    if(not g)DEVABORT("cannot append, not a grid: "+Grid->str());
    std::vector<double> appMesh(_mesh),appWeig(_mesh);
    appMesh.insert(appMesh.end(),g->_mesh.begin(),g->_mesh.end());
    appWeig.insert(appWeig.end(),g->_weights.begin(),g->_weights.end());
    return factory(appMesh,appWeig);
}

const BasisGridQuad* BasisGridQuad::factory(const std::string Definition){
    if(Definition.substr(0,Definition.find(":"))!="GridQuad")
        ABORT("BasisGridQuad definition must start as \"GridQuad:...\", is: "+Definition.substr(0,12)+"....");
    std::vector<std::string> strMesh=tools::splitString(Definition.substr(Definition.find(":")+1),',');
    std::vector<double> meshWeig;
    for(auto s: strMesh)meshWeig.push_back(tools::string_to_double(s));
    return factory(std::vector<double>(meshWeig.begin(),meshWeig.begin()+strMesh.size()/2),
                   std::vector<double>(meshWeig.begin()+strMesh.size()/2,meshWeig.end())  );
}

const BasisGridQuad* BasisGridQuad::factory(const BasisAbstract *Grid){
    BasisGridQuad g(Grid);
    return factory(g.mesh(),g.weights());
}

const BasisGridQuad* BasisGridQuad::factory(std::vector<double> Mesh,std::vector<double> Weights){
    std::vector<double>meshWeig(Mesh);
    meshWeig.insert(meshWeig.end(),Weights.begin(),Weights.end());
    const BasisGridQuad* g=_allBasis[meshWeig];
    if(g==0)g=_allBasis[meshWeig]=new BasisGridQuad(Mesh,Weights);
    return g;
}

std::string BasisGridQuad::strDefinition() const {
    std::string s("GridQuad:");
    for(size_t k=0;k<   _mesh.size();k++)s+=tools::str(   _mesh[k])+",";
    for(size_t k=0;k<_weights.size();k++)s+=tools::str(_weights[k])+",";
    return s.substr(0,s.size()-1); // trim the last ','
}
std::string BasisGridQuad::str(int Level) const{
    std::string s("QuadGrid ");
    s+="["+tools::str(_mesh[0],3)+","+tools::str(_mesh.back(),3)+"] "+tools::str(size())+"["+tools::str(size())+"]";
    return s;
}
