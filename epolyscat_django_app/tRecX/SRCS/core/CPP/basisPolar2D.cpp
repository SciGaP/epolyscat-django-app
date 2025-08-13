// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisPolar2D.h"

#include "axis.h"
#include "index.h"
#include "basisProd.h"
#include "basisAbstract.h"
#include "readInput.h"
#include "printOutput.h"
#include "str.h"
#include "coordinateTrans.h"

using namespace std;

std::vector<double> BasisPolar2D::toCartesian(const std::vector<double>&CoorNdim) const{
    vector<double> vu=CoordinateTrans::fromPolar2d(CoorNdim);
    for(size_t k=0;k<dim();k++)vu[k]+=_origin[k];
    return vu;
}

std::vector<double> BasisPolar2D::fromCartesian(const std::vector<double>&Cartesian) const{
    vector<double> vu(Cartesian);
    for(size_t k=0;k<dim();k++)vu[k]-=_origin[k];
    return CoordinateTrans::toPolar2d(vu);
}
std::vector<std::complex<double> > BasisPolar2D::operator()(std::vector<double> CoorNdim) const{
    return BasisNdim::IndexAt(_idx,CoorNdim);
}

Eigen::MatrixXd BasisPolar2D::jacobianToNdim(const std::vector<double>&CoorNdim) const {
    return Eigen::Map<Eigen::MatrixXd>(CoordinateTrans::jacPolar2d(CoorNdim).data(),2,2);
}

BasisPolar2D::BasisPolar2D(ReadInput & Inp)
{
    _name="Polar2D";
    _ndimCoor="Phi.Rho";

    double rad;
    int mmax,nrad,quad;
    Inp.read("Polar2D","origin",_origin,"0 0","x y, origin of the off-center polar basis");
    Inp.read("Polar2D","radius",rad,"0","radius around origin");
    Inp.read("Polar2D","Mmax",mmax,"0","maximal |m|");
    Inp.read("Polar2D","Nmax",nrad,"10","number of radial functions");
    Inp.read("Polar2D","quad",quad,"10","minimal number of quadrature points");
    if(not Inp.found("Polar2D"))return;
    if(rad==0.)ABORT("must specify radius>0");
    double offRad=sqrt(_origin[0]*_origin[0]+_origin[1]*_origin[1]);
    if(_origin!=vector<double>(2,0.) and rad>0.5*offRad){
        if(rad>=offRad*(1.-1.e-12))
            ABORT("off-center radius must be smaller than distance to origin, is: "+tools::str(rad)+" >= "+tools::str(offRad));
        PrintOutput::warning(
                    Str("off-center basis not well separated from origin, recommended distance > 2*radius, is:")+
                    rad+"vs"+offRad+"\nlarge number of quadrature points may help, is: "+quad);
    }

    // get the main coordinates
    _quadCoor=mainCoor(Inp);
    _idx=idxConstruct(rad,mmax,nrad);

    // get values and first derivatives wrt off-center
    BasisProd off(_idx,quad);
    _valDer=off.valDers();
    _quadGrid=off.quadGrid();
    _quadWeig=off.quadWeig();

    mainQuadValDer();

    test();
}

Index * BasisPolar2D::productIndex(const BasisAbstract *Radial, int Mmax){
    vector<int> marg(2,0);
    marg[1]=2*Mmax;
    BasisSetDef dPhi(2*Mmax+1,0.,2*math::pi,"cosSin",true,true,true,Coordinate("Phi"),marg);

    Index * idx=new Index();
    idx->setAxisName("Phi");
    idx->setBasis(BasisAbstract::factory(dPhi));
    for (size_t l=0;l<idx->basis()->size();l++){
        idx->childAdd(new Index());
        idx->childBack()->setAxisName("Rho");
        idx->childBack()->setBasis(Radial);
        for(int n=0;n<Radial->size();n++)idx->childBack()->leafAdd();
    }
    idx->setFloor(1);
    idx->sizeCompute();
    return idx;

}

std::string BasisPolar2D::selectCoor(const std::string Coor) const{
    vector<vector<string> >allowed;
    allowed.push_back({"X","Y"});
    allowed.push_back({"Phi","Rho"});
    vector<string> coor=tools::splitString(Coor,'.');

    string selCoor;
    if(coor[0]!="Ndim")ABORT("need Ndim as first axis, is: "+Coor);
    for(int k=0;k<allowed.size();k++)
        for(int l=0;l<coor.size();l++)
            if(std::find(allowed[k].begin(),allowed[k].end(),coor[l])!=allowed[k].end())
                selCoor+="."+coor[l];
    if(std::count(selCoor.begin(),selCoor.end(),'.')!=allowed[0].size())
        ABORT("cannot find coordinates to transform to Polar2D: "+Coor);
    return selCoor.substr(1); // remove leading '.'
}

const Index* BasisPolar2D::idxConstruct(double RadMax, int Mmax, int Nrad){
    vector<int> marg(2,0);
    marg[1]=Nrad-1;
    BasisSetDef dR(Nrad,0.,RadMax,"polynomial",true,true,true,Coordinate("Rho"),marg);
    return const_cast<const Index*>(productIndex(BasisAbstract::factory(dR),Mmax));
}
