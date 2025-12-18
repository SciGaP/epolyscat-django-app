// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "mapGauge.h"
#include "discretizationSurface.h"
#include "discretizationGrid.h"
#include "algebra.h"
#include "pulse.h"
#include "coefficientsFloor.h"
#include "basisGrid.h"

#include "operatorIdentity.h"
#include "mpiWrapper.h"

using namespace std;


MapGauge::MapGauge(const DiscretizationSurface *S)
    :OperatorAbstract(S->mapFromParent()->name,S->idx(),S->parent->idx())
{
    //CAUTION: must use new map to grid, old appears to be flawed in some situations (haCC)
    mapToSurface=S->sharedFromParent();
    _construct();
}

MapGauge::MapGauge(const Index* SurfI):OperatorAbstract("GaugeTransform",SurfI,SurfI)
{
    if(SurfI->hierarchy().find("ValDer")==std::string::npos)
        DEVABORT("contructor is for surface discretization, got hierarchy "+SurfI->hierarchy());

    mapToSurface.reset(new OperatorIdentity(SurfI));
    _construct();
}
void MapGauge::_construct(){
    if(not Algebra::isAlgebra("Rg") or Algebra("Rg").val(0.).real()==0)
        DEVABORT("original does not have Gauge radius");

    if(MPIwrapper::Size(MPIwrapper::worldCommunicator())>1)
        DEVABORT("MixedGauge at present only scalar code");
    gridAngular=new DiscretizationGrid(idx(),tools::splitString("Phi.Eta",'.'),
                                       vector<unsigned int>(),vector<vector<double> >(),false);
    double phi(0),eta(0);
    gPhas=new GaugePhase(gridAngular->idx(),phi,eta,rgrid,aDotUnit,phase);

     _cSurf.reset(new Coefficients(mapToSurface->idx()));
     _cGrid.reset(new Coefficients(gridAngular->mapFromParent()->idx()));
}

MapGauge::~MapGauge(){
    delete gridAngular;
    delete gPhas;
}
void MapGauge::apply(std::complex<double> Alfa, const Coefficients &X, std::complex<double> Beta, Coefficients &Y) const {

    mapToSurface->apply(Alfa,X,0.,*_cSurf);
    gridAngular->mapFromParent()->apply(1.,*_cSurf,0.,*_cGrid);
    applyPhase(_cGrid.get(),gPhas);
    gridAngular->mapToParent()->apply(1.,*_cGrid,Beta,Y);
}

void MapGauge::applyPhase(Coefficients *Vec, const GaugePhase *Grid) const {
    // recursion terminates where pPhase is defined
    if(Grid->pPhase!=0){
        *Vec*=*Grid->pPhase;
        // length gauge: angular correction
        if(Grid->pADotUnit!=0)
            for(size_t k=0,l=Vec->size()/2;k<Vec->size()/2;k++,l++)
                Vec->floorData()[l]+=Vec->floorData()[k]*(complex<double>(0.,*Grid->pADotUnit));
    }
    else
        for(size_t k=0;k<Vec->childSize();k++)
            applyPhase(Vec->child(k),Grid->child(k));
}

void MapGauge::update(double Time, const Coefficients* CurrentVec){
    vector<double> a(Pulse::current.vecA(Time)); //watch out: vecA = (Az,Ax,Ay)
    for(size_t k=0;k<rgrid.size();k++){
        aDotUnit[k]=rgrid[k][0]*a[1]+rgrid[k][1]*a[2]+rgrid[k][2]*a[0];
        phase[k]=exp(complex<double>(0.,aDotUnit[k]*rgrid[k][3]));
    }
    if(CurrentVec)return; // shut up about unused variable
}

MapGauge::GaugePhase::GaugePhase(const Index* IdxAng,double Phi, double Eta,
                                 vector<vector<double> > &Unit,deque<double> &ADotUnit,deque<complex<double> > & Phase)
    :pPhase(0),pADotUnit(0){
    double rGauge=Algebra("Rg").val(0.).real();
    if(IdxAng->axisName().find("ValDer")!=string::npos){
        vector<double> coor(4);
        coor[0]=cos(Phi)*sqrt(1-Eta*Eta);
        coor[1]=sin(Phi)*sqrt(1-Eta*Eta);
        coor[2]=Eta;
        coor[3]=min(IdxAng->basis()->grid()->mesh()[0],rGauge);
        size_t k=0;
        while(k<Unit.size() and Unit[k]!=coor)k++;
        if(k==Unit.size()){
            Unit.push_back(coor);
            ADotUnit.push_back(1.);
            Phase.push_back(1.);
        }
        pPhase=&Phase[k];
        if(not IdxAng->hasFloor())ABORT("presently we assume that ValDer is floor level");

        // length-gauge region: angular correction
        pADotUnit=0;
        if(rGauge-IdxAng->basis()->grid()->mesh()[0]>-1.e-10)pADotUnit=&ADotUnit[k];

        return;
    }

    for(size_t k=0;k<IdxAng->childSize();k++){
        if(IdxAng->axisName().find("Phi")==0){
            if(IdxAng->axisName().length()!=3)
                ABORT("illegal phi-coordinate (for now, only single polar coordinats for MapGauge):"+IdxAng->axisName());
            Phi=IdxAng->basis()->grid()->mesh()[k];
        }
        if(IdxAng->axisName()=="Eta")Eta=IdxAng->basis()->grid()->mesh()[k];
        childAdd(new GaugePhase(IdxAng->child(k),Phi,Eta,Unit,ADotUnit,Phase));
    }
}
