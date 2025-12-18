// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "integrateTree.h"
#include "tools.h"

#include "algebra.h"
#include "coefficients.h"
#include "basisIntegrable.h"
#include "index.h"

#include "orthopol.h"

using namespace std;
std::map<std::string,const IntegrationBoundary*> IntegrationBoundary::_list;

IntegrateTree::IntegrateTree(const std::string Definition, const Coefficients &Func, bool SquareFunctions)
    :_definition(Definition),_squareFunctions(SquareFunctions)
{
    // we may not be entering from index top level
    vector<const IntegrationBoundary*>boundaries(Func.idx()->depth(),0);
    for(const Coefficients* f=&Func;f!=0;f=f->descend()){
        boundaries.push_back(IntegrationBoundary::factory(f->idx()->axisName(),Definition));
    }
    while(boundaries.back()==0)boundaries.pop_back();

    vector<int> weightLevels,intLevels;
    for(size_t k=0;k<boundaries.size();k++){
        intLevels.push_back(k);
        if(boundaries[k]!=0)weightLevels.push_back(k);
    }


    _weightIdx.reset(Func.idx()->factor(weightLevels));
    _integralIdx.reset(Func.idx()->factor(weightLevels,true));
    _integralIdx->resetFloor(_integralIdx->heightAboveBottom());

    _weights.reset( new Coefficients(  _weightIdx.get()));
    _integral.reset(new Coefficients(_integralIdx.get()));
    _weights->treeOrderStorage();
    _integral->treeOrderStorage();

    setWeights(&Func,*_weights,{});
}

void IntegrateTree::setWeights(const Coefficients* Integrand, Coefficients & Weights, vector<double> Q){

    if(Weights.isLeaf()){
        // no more integrations, set =1
        Weights.setToConstant(1.);
        return;
    }

    if(not Weights.nodeEquivalent(Integrand)){
        // not integrated over, but further integration levels below
        for(size_t i=0;i<Integrand->childSize();i++)
            setWeights(Integrand->child(i),Weights,Q);
        return;
    }

    // compute the Q-dependent integration rule (qj,wj)
    const IntegrationBoundary* boundary=IntegrationBoundary::factory(Integrand->idx()->axisName(),_definition);
    vector<double> qM,wM;
    // some guess-work about the size of quadrature: exact for multi-variate polynomial behavior
    quadRule(boundary->ranges(Q),boundary->rangeWeights(Q),Weights.size(),qM,wM);

    // evaluate m-level functions for all quadrature points qj
    vector<complex<double > > fM,dum,cqM;
    for(double q: qM)cqM.push_back(q);
    lFunc(Weights)->valDer(cqM,fM,dum);
    // spectrum comes with the wrong functions -> need to introduce modulusSquared-basis
    if(_squareFunctions)
        for(complex<double> & a: fM)a=conj(a)*a;

    // form S[im,im+1,...]=sum[j] fM[im](qj)*wj S[qj;im+1,...] (see Streak.pdf)
    Coefficients weigM(Weights.idx(),0.);
    Q.resize(Q.size()+1);
    for(size_t j=0;j<qM.size();j++){
        Q.back()=qM[j];
        for(size_t i=0;i<Weights.childSize();i++){
            setWeights(Integrand->child(i),*Weights.child(i),Q);
            weigM.child(i)->axpy(fM[j+i*qM.size()]*wM[j],Weights.child(i));
        }
    }
    Weights=weigM;
}

const BasisIntegrable* IntegrateTree::lFunc(const Coefficients &C) const{
    return C.idx()->basis()->integrable();
}

void IntegrateTree::quadRule(const vector<vector<double> >&Ranges, const std::vector<double> &RangeWeights,
                             int Pts, std::vector<double> & Quad, std::vector<double> &Weig)
{
    if(Ranges.size()>RangeWeights.size())
        DEVABORT("there are fewer interval weights than intervals");
    // several disconnected intervals
    Quad.clear();
    Weig.clear();

    for(size_t k=0;k<Ranges.size();k++){
        vector<double>q,w;
        quadRule(Ranges[k],RangeWeights[k],Pts,q,w);
        Quad.insert(Quad.end(),q.begin(),q.end());
        Weig.insert(Weig.end(),w.begin(),w.end());
    }
}

// for now: Legendre quadratur
void IntegrateTree::quadRule(const vector<double> &Range, double RangeWeight, int Pts, std::vector<double> & Quad, std::vector<double> &Weig){
    OrthogonalLegendre leg;
    leg.quadrature(Pts,Quad,Weig);
    for(size_t k=0;k<Quad.size();k++){
        Weig[k]*=(Range[1]-Range[0])/2.*RangeWeight;
        Quad[k]=Range[0]+(1+Quad[k])*(Range[1]-Range[0])/2.;
    }
}

Coefficients IntegrateTree::integrate(const Coefficients* Integrand){
    Coefficients Integral(_integralIdx.get());
    Integral.setToZero();
    integrate(Integral,Integrand,_weights.get());
    return Integral;
}


void IntegrateTree::integrate(Coefficients& Integral, const Coefficients *Integrand, const Coefficients* Weight)
{
    if(Weight->isLeaf())
    { // no more integrations
        Integral.axpy(Weight->anyData()[0],Integrand);
    }
    else if(not Integrand->nodeEquivalent(Weight)){
        // do not integrate on this level
        for(size_t k=0;k<Integrand->childSize();k++)
            integrate(*Integral.child(k),Integrand->child(k),Weight);
    }
    else {
        for(size_t k=0;k<Integrand->childSize();k++){
            integrate(Integral,Integrand->child(k),Weight->child(k));
        }
    }
}

const IntegrationBoundary* IntegrationBoundary::factory(string Axis, string Definition){
    string axDef=Axis+":"+Definition;
    if(not _list.count(axDef)){
        string kind=Definition.substr(0,Definition.find("["));
        vector<string> part=tools::splitString(tools::stringInBetween(Definition,"[","]"),',');
        vector<double> angl;
        if(part.front()!=Definition)
            for(string p: part)angl.push_back(Algebra::realConstant(p));

        if(kind=="cone"){
            if(angl.size()!=3)
                ABORT("need definition of cone axis and opening angle as cone[phi,theta,gamma], have: "+Definition);
            if(angl[2]<0.)ABORT("need non-negative cone angle, have: "+Definition);
            if(Axis.find("Phi")==0)
                _list[axDef]=new IntegrationConePhi(angl[0],angl[1],angl[2]);
            else if(Axis.find("Eta")==0)
                _list[axDef]=new IntegrationConeEta(angl[0],angl[1],angl[2]);
            else
                _list[axDef]=0;
        }

        else if(kind=="zone"){
            if(angl.size()!=2)
                ABORT("need definition as zone[thetaMin,thetaMax], have: "+Definition);
            if(Axis.find("Eta")==0){
                _list[axDef]=new IntegrationZoneEta(angl[0],angl[1]);
            }
            else
                _list[axDef]=0;
        }

        else
            ABORT("no DomainBoundary defined for type "+kind+" of "+Definition);
    }
    return _list[axDef];
}
IntegrationZoneEta::IntegrationZoneEta(double LowerTheta,double UpperTheta){
    double width=UpperTheta-LowerTheta;
    if(width>math::pi+1.e-12)ABORT(Str("zone width must be <= pi, is: ")+(width/math::pi)+"pi");
    if(LowerTheta>=0. and UpperTheta<=(math::pi+1.e-12))
        _ranges={{cos(LowerTheta),cos(UpperTheta)}};
    else if (LowerTheta<0.)
        _ranges={{cos(LowerTheta),1.},{cos(UpperTheta),1.}};
    else if (UpperTheta>(math::pi+1.e-12))
        _ranges={{-1.,cos(LowerTheta)},{-1.,cos(UpperTheta)}};
    else
        DEVABORT("this is not supposed to happen: "+tools::str(LowerTheta)+" "+tools::str(UpperTheta));
}

IntegrationConePhi::IntegrationConePhi(double Phi, double Theta, double Gamma){
    if(Theta<Gamma and math::pi-Theta<Gamma){
        // both poles in range - probe complement separately
        IntegrationConePhi compPhi(Phi+math::pi,math::pi-Theta,math::pi-Gamma);
        _ranges=compPhi.ranges({});
        _ranges.push_back({compPhi.ranges({})[0][1],compPhi.ranges({})[0][0]+2*math::pi});
    }
    else if(Theta<Gamma or math::pi-Theta<Gamma or Gamma>math::pi/2){
        // one of the poles is in the cone, whole range
        _ranges.push_back({0.,2*math::pi});
    }
    else {
        double eta=cos(Theta)/cos(Gamma);
        double arg=(cos(Gamma)-cos(Theta)*eta)/(sin(Theta)*sqrt(1-eta*eta));
        if(std::abs(arg)>1.){
            if(abs(arg)<1.e-12)arg=(arg>0) ? 1. : -1.;
            DEVABORT("incorrect argument");
        }
        _ranges.push_back({Phi-acos(arg),Phi+acos(arg)});
    }
}


// this should go into namespace tools::
void snapValue(double &Value, std::vector<double> SnapGrid,double Eps=1.e-12){
    double gmin=DBL_MAX,gmax=-DBL_MAX;
    for(double g: SnapGrid){
        gmin=min(g,gmin);
        gmax=max(g,gmin);
    }
    for(double g: SnapGrid)
        if(std::abs(g-Value)<Eps*(gmax-gmin))Value=g;
}

IntegrationConeEta::IntegrationConeEta(double Phi, double Theta, double Gamma)
    :_phiComp(0),_etaComp(0),theta(Theta),gamma(Gamma)
{
    // move exact boundary points if small differencies
    snapValue(theta,{0.,math::pi});
    snapValue(gamma,{0.,math::pi});

    if(theta>math::pi and theta<math::pi+1.e12)theta=math::pi;

    if(gamma <0. or math::pi<Gamma)ABORT(Str("choose cone gamC in [0,pi], is: ")+(Gamma/math::pi)+"pi");
    if(theta <0. or theta>math::pi)ABORT(Str("choose polar axis theta in [0,pi], is: ")+(Theta/math::pi)+"pi");
    cz=cos(Theta);
    cx=sin(Theta)*cos(Phi);
    cy=sin(Theta)*sin(Phi);
    pPart=-2*cos(Gamma)*cz;
    qPart=std::pow(cos(Gamma),2);

    if(theta<gamma and math::pi-theta<gamma){
        // both poles in domain - boundaries of complement
        std::string comp=Str("cone,","")+(Phi+math::pi)+","+(math::pi-Theta)+","+(math::pi-Gamma);
        _phiComp=IntegrationBoundary::factory("Phi",comp);
        _etaComp=IntegrationBoundary::factory("Eta",comp);
    }
}

const vector<vector<double> > IntegrationConeEta::ranges(const std::vector<double> &Q)const{
    if(Q.size()!=1)DEVABORT(Str("need exactly 1 paramter, found")+Q.size());
    if(_phiComp!=0){ // both poles in domain

        // not in phi-range of complement
        if(Q[0]<=_phiComp->ranges({})[0][0] or _phiComp->ranges({})[0][1]<=Q[0])
            return {{-1.,1.}};

        // two intervals: from south pole to complement, from complement to north pole
        vector<double> compRange=_etaComp->ranges(Q)[0];
        return {{-1.,compRange[0]},{compRange[1],1.}};
    }

    double aSqu=std::pow(cx*cos(Q[0])+cy*sin(Q[0]),2);
    double R=aSqu+cz*cz;
    double p2=pPart/(2*R);
    double q=(qPart-aSqu)/R;
    double root=p2*p2-q;
    if(root<0){
        if(abs(root)>1.e-12)
            DEVABORT(Str("error: X,root,x,y,z,z*cos(Gamma)=")+Q+root+cx+cy+cz+(-qPart*0.5));
        root=0.;
    }
    if(theta>=gamma and math::pi-theta>=gamma)return {{-p2-sqrt(root),-p2+sqrt(root)}}; // neither pole in domain
    if(theta <gamma and math::pi-theta>=gamma)return {{-p2-sqrt(root), 1.}}; // only north pole in domain
    if(theta>=gamma and math::pi-theta <gamma)return {{-1.,-p2-sqrt(root)}}; // only south pole in domain
    DEVABORT("failure - check code");

}

