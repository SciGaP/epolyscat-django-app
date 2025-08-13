// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisPolarOff.h"

#include "readInput.h"
#include "printOutput.h"
#include "axis.h"
#include "index.h"
#include "basisAbstract.h"
#include "coordinateTrans.h"
#include "basisProd.h"
#include "qtEigenDense.h"
#include "str.h"
//#include "basisMat.h"
#include "axis.h"
#include "tools.h"

using namespace std;

std::vector<std::complex<double> > BasisPolarOff::operator()(std::vector<double> CoorOffPolar) const{
    return IndexAt(_idx,CoorOffPolar);
}

bool BasisPolarOff::hasOverlap(const BasisAbstract &Other) const{
    const BasisPolarOff * o=dynamic_cast<const BasisPolarOff*>(&Other);
    if(not o)DEVABORT("cannot determine overlap between basis functions "+name()+" and "+Other.name());
    double dist(0.);
    for(size_t i=0;i<3;i++)dist+=std::pow(_origin[i]-o->_origin[i],2);
    double rThis=tools::string_to_double(tools::splitString(strDefinition(),',')[3]);
    double rOther=tools::string_to_double(tools::splitString(strDefinition(),',')[3]);
    return std::pow(rThis+rOther,2)>dist;
}

std::vector<double> BasisPolarOff::toCartesian(const std::vector<double>&CoorNdim) const{
    if(_ndimCoor!="Phi.Eta.Rn")DEVABORT("only implemented for Phi.Eta.Rn");
    vector<double> vu=CoordinateTrans::fromPolar3d(CoorNdim);
    for(size_t k=0;k<dim();k++)vu[k]+=_origin[k];
    return vu;
}

std::vector<double> BasisPolarOff::fromCartesian(const std::vector<double>&XYZ) const{
    vector<double> xyz(XYZ);
    for(size_t k=0;k<3;k++)xyz[k]-=_origin[k];
    return CoordinateTrans::toPolar3d(xyz);
}

Eigen::MatrixXd BasisPolarOff::jacobianToNdim(const std::vector<double>&CoorNdim) const {
    return Eigen::Map<Eigen::MatrixXd>(CoordinateTrans::jacPolar3d(CoorNdim).data(),3,3);
}

BasisPolarOff::BasisPolarOff(ReadInput &Inp,std::string Category,int Line)
{
    _name="PolarOffCenter"+tools::str(Line,2,'0');
    _ndimCoor="Phi.Eta.Rn";
    double rad,rmin;
    int lmax,mmax,nrad,quad;
    std::string sOrig;
    Inp.texdocuCategoryAdd(Category,"origin,radius,Nmax,Mmax,Lmax,quad",
                           R"tex(
                           Off-center basis $\Phi_i(\phi,\eta,r)$ in spherical coordinates
                           consisting of spherical harmonics and polynomials in $r$.
                           Useful for capturing, e.g., off-center atoms in a molecule.
                           )tex","510");
    Inp.read(Category,"origin",sOrig,"[0,0,0]","[x,y,z], origin of the off-center polar basis",Line);
    Inp.read(Category,"radius",rad,"0","radius around origin",Line)
            .texdocu(R"tex(
                     Basis has support on $r=[\text{rmin},\text{radius})$, =0 at boundaries
                     )tex");
    Inp.read(Category,"rmin",rmin,"0","radial support is in [rmin,radius]",Line)
            .texdocu(R"tex(
                     Basis has support on $r=[\text{rmin},\text{radius})$, =0 at radius.
                     )tex");
    Inp.read(Category,"Lmax",lmax,"2","maximal angular momentum",Line)
            .texdocu(R"tex(
                     Use spherical harmonics $Y_{lm}$ up to $l\leq$\lcode{Lmax}.
                     )tex");
    Inp.read(Category,"Mmax",mmax,tools::str(lmax),"maximal |m|",Line)
            .texdocu(R"tex(
                     Use spherical harmonics $Y_{lm}$ up to $|m|\leq$\lcode{Mmax}.
                     )tex");
    Inp.read(Category,"Nmax",nrad,"10","number of radial functions",Line)
            .texdocu(R"tex(
                     Basis is $Y_{lm}(\phi,\eta)P_n(r)$ with polynomials of degree $n<$\lcode{Nmax}.
                     )tex");
    Inp.read(Category,"quad",quad,"10","minimal number of quadrature points",Line)
            .texdocu(R"tex(
                     A brute-force quadrature will be performed with the main discritization using
                     high order quadratures. This sets a lower bound for the number of quadrature points.
                     )tex");
    if(not Inp.found(Category) or (rad==0. and Inp.endCategory(Category,Line)))return; // no radius specified
    if(std::max(0.,rad-rmin)<=rad*0.1)PrintOutput::warning(Sstr+"very thin shell in Ndim polar basis: [rmin,rmax]=["+rmin+","+rad+"]");
    if(rad<=0.)ABORT("must specify radius>0"); // non-positive radius


    if(sOrig.find(",")==std::string::npos
            or sOrig.find("[")==std::string::npos
            or sOrig.find("]")==std::string::npos
            )ABORT("input format origin="+sOrig+"is obsolete, use origin=[x,y,z] instead");
    _origin.clear();
    _strDefinition="PolarOffCenter:"+sOrig+","+tools::str(rad)+","+tools::str(lmax)+","+tools::str(mmax)+","+tools::str(nrad)+","+tools::str(quad);
    for(auto s: tools::splitString(tools::stringInBetween(sOrig,"[","]"),','))
        _origin.push_back(tools::string_to_double(s));

    double offRad=sqrt(_origin[0]*_origin[0]+_origin[1]*_origin[1]+_origin[2]*_origin[2]);
    if(_origin!=vector<double>(3,0.) and rad>0.5*offRad){
        if(rad>=offRad*(1.-1.e-12))
            ABORT("off-center radius must be smaller than distance to origin, is: "+tools::str(rad)+" >= "+tools::str(offRad));
        PrintOutput::warning(
                    Str("off-center basis not well separated from origin, recommended distance > 2*radius, is:")+
                    rad+"vs"+offRad+"\nlarge number of quadrature points may help, is: "+quad);
    }

    // get the reference coordinates
    _quadCoor=mainCoor(Inp);
    _idx=idxConstruct(rmin,rad,lmax,mmax,nrad);

    // get native quadrature points and values and first derivatives there
    BasisProd off(_idx,quad);
    _valDer=off.valDers();
    _quadGrid=off.quadGrid();
    _quadWeig=off.quadWeig();

    // transform to main coordinates and cooresponding values and derivatives
    mainQuadValDer();


    test();
}

const Index* BasisPolarOff::idxConstruct(double RadMin, double RadMax, int Lmax, int Mmax, int Nrad){
    vector<int> marg(2,0);
    marg[1]=Nrad-1;
    BasisSetDef dR(Nrad,RadMin,RadMax-RadMin,"polynomial",true,true,true,Coordinate("Rn"),marg);
    const BasisAbstract *rBas=BasisAbstract::factory(dR);
    return const_cast<const Index*>(productIndex(rBas,Lmax,Mmax));
}

Index * BasisPolarOff::productIndex(const BasisAbstract *Radial, int Lmax, int Mmax){
    vector<int> marg(2,0);
    marg[1]=2*Mmax;
    BasisSetDef dPhi(2*Mmax+1,0.,2*math::pi,"cosSin",true,true,true,Coordinate("Phi"),marg);

    Index * idx=new Index();
    idx->setAxisName("Phi");
    idx->setBasis(BasisAbstract::factory(dPhi));
    for(size_t k=0;k<idx->basis()->size();k++){
        marg[1]=Lmax-abs(idx->basis()->physical(k));
        BasisSetDef dEta(Lmax+1-abs(idx->basis()->physical(k)),-1.,2.,"assocLegendre",
                         true,true,true,Coordinate::fromString("Eta"),marg,
                         ComplexScaling(),false,vector<double>(1,idx->basis()->physical(k)));
        idx->childAdd(new Index());
        idx->childBack()->setAxisName("Eta");
        idx->childBack()->setBasis(BasisAbstract::factory(dEta));
        for (size_t l=0;l<idx->childBack()->basis()->size();l++){
            idx->childBack()->childAdd(new Index());
            idx->childBack()->childBack()->setAxisName("Rn");
            idx->childBack()->childBack()->setBasis(Radial);
            for(size_t n=0;n<Radial->size();n++)idx->childBack()->childBack()->leafAdd();
        }
    }
    idx->setFloor(2);
    idx->sizeCompute();
    return idx;

}
