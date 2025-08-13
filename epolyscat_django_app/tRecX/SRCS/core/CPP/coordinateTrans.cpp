// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coordinateTrans.h"
#include "qtEigenDense.h"
#include <complex>
#include "constants.h"
#include "str.h"
#include "abort.h"


using namespace std;

// things below should be a single routine
CoordinateMap CoordinateTrans::toCartesian(string Coor){
    if      (Coor=="X.Y")       return identity;
    else if (Coor=="Phi.Rho")    return fromPolar2d;
    else if (Coor=="X.Y.Z")      return identity;
    else if (Coor=="Phi.Eta.Rn") return fromPolar3d;
    else if (Coor=="Phi.Eta.kRn") return fromPolar3d;
    else if (Coor=="Phi.Rho.Z") return fromCylinder3d;
    else return undefined;
}

CoordinateMap CoordinateTrans::fromCartesian(string Coor){
    if      (Coor=="X.Y")        return identity;
    else if (Coor=="Phi.Rho")    return toPolar2d;
    else if (Coor=="X.Y.Z")      return identity;
    else if (Coor=="Phi.Eta.Rn") return toPolar3d;
    else if (Coor=="Phi.Eta.kRn") return toPolar3d;
    else if (Coor=="Phi.Rho.Z") return toCylinder3d;
    else return undefined;
}
IntegrationFactor CoordinateTrans::integrationFactor(string Coor){
    if      (Coor=="X.Y")        return constOne;
    else if (Coor=="Phi.Rho")    return constOne;
    else if (Coor=="X.Y.Z")      return constOne;
    else if (Coor=="Phi.Eta.Rn") return radius;
    else ABORT("undefined transformation from Cartesian to "+Coor);
}
NablaFactor CoordinateTrans::nablaFactor(string Coor){
    if      (Coor=="X.Y")        return nablaOne;
    else if (Coor=="Phi.Rho")    return nablaOne;
    else if (Coor=="X.Y.Z")      return nablaOne;
    else if (Coor=="Phi.Eta.Rn") return nablaRadius;
    else ABORT("undefined transformation from Cartesian to "+Coor);
}
JacobianMap CoordinateTrans::jacCartesian(string Coor){
    if      (Coor=="X.Y")        return jacIdentity;
    else if (Coor=="Phi.Rho")    return jacPolar2d;
    else if (Coor=="X.Y.Z")      return jacIdentity;
    else if (Coor=="Phi.Eta.Rn") return jacPolar3d;
    else if (Coor=="Phi.Eta.kRn") return jacPolar3d;
    else if (Coor=="Phi.Rho.Z") return jacCylinder3d;
    else ABORT("undefined jacobian to Cartesian from "+Coor);
}

std::vector<double> CoordinateTrans::jacIdentity(const std::vector<double> & In){
    if(In.size()==2)return std::vector<double>({1.,0.,0.,1.});
    if(In.size()==3)return std::vector<double>({1.,0.,0.,0.,1.,0.,0.,0.,1.});
    ABORT(Str("no jacIdentity defined for In.size()==")+In.size());
    return {};
}

vector<double> CoordinateTrans::shiftPolar(std::vector<double> PhiEtaR, std::vector<double> CartesianShift){
    vector<double> xyz(fromPolar3d(PhiEtaR));
    xyz[0]+=CartesianShift[0];
    xyz[1]+=CartesianShift[1];
    xyz[2]+=CartesianShift[2];
    return toPolar3d(xyz);
}
vector<complex<double> > CoordinateTrans::shiftPolar(std::vector<std::complex<double> > PhiEtaR, std::vector<double> CartesianShift){
    vector<double> realCoor(3);
    for(int k=0;k<3;k++)realCoor[k]=PhiEtaR[k].real();
    realCoor=shiftPolar(realCoor,CartesianShift);
    vector<complex<double> > cplxCoor(PhiEtaR);
    for(int k=0;k<3;k++)cplxCoor[k].real(realCoor[k]);
    return cplxCoor;
}
vector<double> CoordinateTrans::fromPolar3d(const std::vector<double> &In){
    double sinPhi=sin(In[0]);
    double cosPhi=cos(In[0]);
    double cosThe=In[1];
    double sinThe=sqrt(1.-In[1]*In[1]);

    vector<double> xyz;
    xyz.push_back(In[2]*sinThe*cosPhi);
    xyz.push_back(In[2]*sinThe*sinPhi);
    xyz.push_back(In[2]*cosThe);
    return xyz;
}

vector<double> CoordinateTrans::toCylinder3d(const std::vector<double> &In){

    vector<double> phiRhoZ;
    double rho=sqrt(In[0]*In[0]+In[1]*In[1]);
    if(rho==0.)return {0.,0.,In[2]};

    // note: the two-argument atan takes care of quadrant
    phiRhoZ.push_back(atan2(In[1],In[0]));
    if(phiRhoZ[0]<0)phiRhoZ[0]+=2*math::pi;
    phiRhoZ.push_back(rho);
    phiRhoZ.push_back(In[2]);
    return phiRhoZ;
}
vector<double> CoordinateTrans::fromCylinder3d(const std::vector<double> &In){
    double sinPhi=sin(In[0]);
    double cosPhi=cos(In[0]);

    vector<double> xyz;
    xyz.push_back(In[1]*cosPhi);
    xyz.push_back(In[1]*sinPhi);
    xyz.push_back(In[2]);
    return xyz;
}

vector<double> CoordinateTrans::jacCylinder3d(const std::vector<double> &PhiRhoZ){
    double sinPhi=sin(PhiRhoZ[0]);
    double cosPhi=cos(PhiRhoZ[0]);

    vector<double> res;
    res.push_back(-sinPhi*PhiRhoZ[1]);
    res.push_back(cosPhi);
    res.push_back(0.);

    res.push_back( cosPhi*PhiRhoZ[1]);
    res.push_back(sinPhi);
    res.push_back(0.);

    res.push_back(0.);
    res.push_back(0.);
    res.push_back(1.);
    return res;
}


vector<double> CoordinateTrans::toPolar3d(const std::vector<double> &In){

    vector<double> phiEtaR;
    double r=sqrt(In[0]*In[0]+In[1]*In[1]+In[2]*In[2]);
    if(r==0.)return vector<double>(3,0.);

    // note: the two-argument atan takes care of quadrant
    phiEtaR.push_back(atan2(In[1],In[0]));
    if(phiEtaR[0]<0)phiEtaR[0]+=2*math::pi;
    phiEtaR.push_back(In[2]/r);
    phiEtaR.push_back(r);
    return phiEtaR;
}

vector<double> CoordinateTrans::fromPolar2d(const std::vector<double> &In){

    vector<double> xy;
    xy.push_back(In[1]*cos(In[0]));
    xy.push_back(In[1]*sin(In[0]));
    return xy;
}
vector<double> CoordinateTrans::toPolar2d(const std::vector<double> &In){

    double rho=sqrt(In[0]*In[0]+In[1]*In[1]);
    if(rho==0.)return vector<double>(2,0.);

    vector<double> phiRho;
    // note: the two-argument atan takes care of quadrant
    phiRho.push_back(atan2(In[1],In[0]));
    if(phiRho[0]<0)phiRho[0]+=2*math::pi;
    phiRho.push_back(rho);
    return phiRho;
}

vector<double> CoordinateTrans::jacPolar3d(const std::vector<double> &PhiEtaR){
    double sinPhi=sin(PhiEtaR[0]);
    double cosPhi=cos(PhiEtaR[0]);
    double cosThe=PhiEtaR[1];
    double sinThe=sqrt(1.-PhiEtaR[1]*PhiEtaR[1]);

    vector<double> res;
    res.push_back(-sinPhi*sinThe*PhiEtaR[2]);
    res.push_back(-cosPhi*cosThe/sinThe*PhiEtaR[2]);
    res.push_back(cosPhi*sinThe);

    res.push_back( cosPhi*sinThe*PhiEtaR[2]);
    res.push_back(-sinPhi*cosThe/sinThe*PhiEtaR[2]);
    res.push_back(sinPhi*sinThe);

    res.push_back(0.);
    res.push_back(PhiEtaR[2]);
    res.push_back(cosThe);

    return res;
}

vector<double> CoordinateTrans::jacPolar2d(const std::vector<double> &PhiRho){
    double sinPhi=sin(PhiRho[0]);
    double cosPhi=cos(PhiRho[0]);

    vector<double> res;
    res.push_back(-sinPhi*PhiRho[1]);
    res.push_back(cosPhi);

    res.push_back( cosPhi*PhiRho[1]);
    res.push_back(sinPhi);
    return res;
}
