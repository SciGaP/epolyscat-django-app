// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "fieldTailor.h"
#include "pulse.h"
#include "orthopol.h"
#include "readInputRange.h"
#include "printOutput.h"
#include "asciiFile.h"
#include "constants.h"
#include "algebra.h"

using namespace std;


std::vector<std::string> PhaseSpace::compNames;
void PhaseSpace::setup(){
    if(compNames.size()>0)return;
    compNames.push_back("Rz");
    compNames.push_back("Rx");
    compNames.push_back("Ry");
    compNames.push_back("Pz");
    compNames.push_back("Px");
    compNames.push_back("Py");
}

bool PhaseSpace::inGoing(const std::vector<double> & Poff) const {
    double s=data()[0]*(data()[3]-Poff[0])+data()[1]*(data()[4]-Poff[1])+data()[2]*(data()[5]-Poff[2]);
    return s<0.;
}


FieldTailor::FieldTailor(ReadInput & Inp)
{
    bool plotPulse;
    Inp.read("FieldTailor","printPulse",plotPulse,"true","write laser to a file");
    Pulse::read(Inp,plotPulse);
    parRange=ReadInputRange(Inp);

    tMax=Units::convert(2.,"OptCyc");
    tMin=Units::convert(1.,"OptCyc");
    Inp.read("FieldTailor","tMin",tMin,tools::str(tMin,3),"minimal time-propagation");
    Inp.read("FieldTailor","tMax",tMax,tools::str(tMax,3),"maximal time-propagation");
    Inp.read("FieldTailor","rescatterRadius",rsqApproach,Algebra::Infinity,"radius considered as re-scattering");
    if(rsqApproach<sqrt(DBL_MAX/2.))rsqApproach*=rsqApproach*100;

    Inp.read("FieldTailor","kind",kind,"recollide","possible criteria: recollide");
    Inp.read("FieldTailor","hamiltonian",ham,"<<FreeMotionInDipole>>","Hamiltonian function for classical motion");

    Inp.read("FieldTailor","printComp",writeComp,"Rx Ry Rz Px Py Pz","print these phase-space components");

    if(ham=="<<FreeMotionInDipole>>")dyn=new FreeMotionInDipole(&Pulse::current,8);
    else ABORT("undefined hamiltonian: "+ham);

    outDir=Inp.outputTopDir();
    maxStep=Units::convert(1.,"OptCyc")/10.;
    tol=1.e-8;


    Inp.read("FieldTailor","optimize",optimize,"true","optimize after grid search",1,"optimize");

    Inp.read("FieldTailor","recollisionPosition",positionTarget,"0 0 0","where to direct trajectory");
    Inp.read("FieldTailor","recollisionEnergy",momentumSquTarget,"0","where to direct trajectory");
    momentumSquTarget*=2.;
    Inp.read("FieldTailor","recollisionAngle",angleTarget,"0","desired angle (in degrees) for recollision");
    angleTarget*=math::pi/180.;

    Inp.read("FieldTailor","sigmaEnergy",sigmaMomentumSqu,"1.","width for momentum square distribution");
    sigmaMomentumSqu*=2.;
    Inp.read("FieldTailor","sigmaAngle",sigmaAngle,"3","width for angular deviation (degrees)");
    sigmaAngle*=math::pi/180.;

    if(not Inp.found("FieldTailor","recollisionAngle"))sigmaAngle=0.;
    if(not Inp.found("FieldTailor","sigmaEnergy"))sigmaMomentumSqu=0.;

    // limits and initial region for optimization
    lowLimit=dlib::matrix<double>(parRange.size(),1);
    uppLimit=dlib::matrix<double>(parRange.size(),1);
    rhoBegin=DBL_MAX;
    for(int k=0;k<parRange.size();k++){
        lowLimit(k)=parRange.lowVal(k);
        uppLimit(k)=parRange.upVal(k);
        rhoBegin=min(rhoBegin,(uppLimit(k)-lowLimit(k)));
    }
    rhoBegin*=0.3;
}

void FieldTailor::print() const {
    Pulse::current.output();
    parRange.print();
    PrintOutput::newLine();
    PrintOutput::lineItem("min",tMin);
    PrintOutput::lineItem("max time",tMax);
    PrintOutput::newLine();
    PrintOutput::lineItem("recoll. radius",sqrt(rsqApproach));
    PrintOutput::lineItem("angle",angleTarget*180./math::pi);

    PrintOutput::paragraph();
    PrintOutput::paragraph();
}

void FieldTailor::run(){

    AsciiFile f(outDir+"pars");
    vector<string>comm;
    comm.push_back("");
    for(int k=0;k<parRange.size();k++)
        comm.back()+=string(" "+parRange.name(k)+"-"+tools::str(parRange.line(k))+"            ").substr(0,13);
    comm.back()+=" tstart      tstop       F(t0)       rmin        E_coll      angle_coll";
    f.writeComments(comm);

    vector<double> par;
    int total=0;
    cout<<"running ";

    while (parRange.next(par))
    {
        total++;
        for(int k=0;k<par.size();k++)
            Pulse::current.updateParameter(parRange.line(k)-1,parRange.name(k),par[k]);
        double tstart=startTime();
        Pulse::current.alignFieldZ(tstart);

        traj.push_back(runSingle(tstart,PhaseSpace()));

        // keep only within time limit
        if(timeExceeded(traj.back())){
            traj.pop_back();
        }
        else {
            if(traj.back().uVec.back().normSqu()>rsqApproach)
                traj.pop_back();
            else
            {
                vector<double>opt(par);
                if(optimize){
                    pulseOptimize(opt);
                }
                // rerun trajectory with new parameters and small step
                double defaultStep=maxStep;
                maxStep=Units::convert(0.01,"OptCyc");
                traj.back()=runSingle(startTime(),PhaseSpace());
                maxStep=defaultStep;

                if(optimize and
                        (abs(angle(traj.back().time.back(),traj.back().uVec.back().p())-angleTarget)>sigmaAngle
                        or traj.back().uVec.back().q().normSqu()>rsqApproach))
                {
                    cout<<".";
                    traj.pop_back();
                }
                else {
                    writeData(f,opt,traj.back());
                    if(traj.size()%50==0)cout<<"\n"<<traj.size();
                    cout<<"+";
                }
            }
        }
    }
    cout<<endl;

    PrintOutput::message(tools::str(traj.size())+" out of "+tools::str(total)+" trajectories meet criterion");
    PrintOutput::message("pulse parameters and trajectories on file "+f.name());
}


void FieldTailor::writeData(const AsciiFile &f,std::vector<double> Par, const Trajectory &Traj){

    vector<double> par(Par);
    par.push_back(Traj.time[0]);
    par.push_back(Traj.time.back());

    // field at start time
    par.push_back(Pulse::current.vecF(Traj.time[0])[0]);

    // add recollision distance (au), energy (au), and angle (degrees)
    VectorReal q=Traj.uVec.back().q();
    VectorReal p=Traj.uVec.back().p();
    par.push_back(sqrt(q.normSqu()));
    par.push_back(0.5*(p-Pulse::current.vecA(Traj.time.back())).normSqu());
    par.push_back(angle(Traj.time.back(),p)*180./math::pi);
    VectorReal a(Pulse::current.vecA(Traj.time.back()));
    par.push_back(-a[0]);
    par.push_back(-a[1]);

    for(int l=0;l<parRange.size();l++){
        if(parRange.name(l).find("(")!=string::npos){
            par[l]=Units::convert(par[l],"DEFAULT_SYSTEM",tools::stringInBetween(parRange.name(l),"(",")"));
        }
    }
    f.writeRow(par,4);
}


double FieldTailor::startTime() const{
    // primitive start time - peak-field
    double t=Pulse::current.beginPrint();
    double dt=(Pulse::current.endPrint()-t)/1000.;
    double fmax=Pulse::current.Fabs(t).real(),tStart=t;
    for(int k=0;k<1000;k++,t+=dt){
        if(fmax<Pulse::current.Fabs(t).real()){
            fmax=Pulse::current.Fabs(t).real();
            tStart=t;
        }
    }
    return tStart;
}

void FieldTailor::zeroIn(Trajectory &Traj) const {
    // find closest approach
    int last=Traj.uVec.size()-1;

    while((Traj.uVec[last]-Traj.uVec[last-1]).normSqu()>0.001){

        double tHalf=0.5*(Traj.time[last]-Traj.time[last-1]);
        PhaseSpace uVec=Traj.uVec[last-1];
        dyn->step(uVec,Traj.time[last-1],tHalf);
        if(uVec.inGoing(Pulse::current.vecA(Traj.time[last-1]+tHalf,false))){
            Traj.uVec[last-1] =uVec;
            Traj.time[last-1]+=tHalf;
        } else {
            Traj.uVec[last] =uVec;
            Traj.time[last]-=tHalf;
        }
    }
    return;
}

void FieldTailor::pulseOptimize(std::vector<double>  &pars) const{
    column_vector start(pars.size());
    // move to within limits (violation may occur due to roundoff)
    for(int k=0;k<pars.size();k++)
        start(k)=max(lowLimit(k),min(uppLimit(k),pars[k]));
    find_min_bobyqa(*this,start,
                    9,    // number of interpolation points
                    lowLimit,uppLimit,
                    rhoBegin,    // initial trust region radius
                    rhoBegin*1.e-2,  // stopping trust region radius
                    10000    // max number of objective function evaluations
                    );
    for(int k=0;k<pars.size();k++)pars[k]=start(k);
}

void FieldTailor::writeTrajectories() const {

    vector<string> wCom=tools::splitString(writeComp,' ');
    vector<int> wNum;
    for(int k=0;k<wCom.size();k++){
        wNum.push_back(0);
        while(wNum.back()<PhaseSpace::compNames.size()
              and PhaseSpace::compNames[wNum.back()]!=tools::cropString(wCom[k])
              )wNum.back()++;
        if(wNum.back()==PhaseSpace::compNames.size())
            ABORT("no phase space component \""+wCom[k]
                  +"\", admissible: "+tools::str(PhaseSpace::compNames));
    }

    vector<vector<double> > cols;
    for(int k=0;k<traj.size();k++){
        cols.push_back(vector<double>());
        for (int l=0;l<traj[k].time.size();l++){
            cols.back().push_back(traj[k].time[l]);
        }
        for (int m=0;m<wNum.size();m++){
            cols.push_back(vector<double>());
            for (int l=0;l<traj[k].time.size();l++)
                cols.back().push_back(traj[k].uVec[l][wNum[m]]);
        }
    }

    AsciiFile f(outDir+"traj",INT_MAX,", ",false);
    vector<string> comm(1,"time "+writeComp);
    f.writeComments(comm);
    f.writeCols(cols,5);
    PrintOutput::message("trajectory components " +tools::str(writeComp)+" on file "+f.name());

}

FieldTailor::Trajectory FieldTailor::runSingle(double T0, const PhaseSpace &UStart) const {
    Trajectory tra;
    tra.uVec.push_back(UStart);
    tra.time.push_back(T0);
    int iStop=0;
    double step=maxStep*0.1;
    PhaseSpace errVec;
    while(not stopCriterion(tra,iStop)){
        dyn->stepError(tra.uVec,tra.time.back(),step,errVec);
        tra.time.push_back(tra.time.back()+step*0.5);
        tra.time.push_back(tra.time.back()+step*0.5);
        if(not dyn->acceptStep(tol,step,step,errVec.norm())){
            tra.uVec.resize(tra.uVec.size()-2);
            tra.time.resize(tra.time.size()-2);
        }
        if(step<1.e-11)ABORT("step underflow");
        step=min(maxStep,step);
        if(tra.time.size()<100)step=min(step,0.1);
    }
    tra.uVec.resize(iStop+1);
    tra.time.resize(iStop+1);

    if(not timeExceeded(tra))zeroIn(tra);
    return tra;
}

bool FieldTailor::isRecolliding(const Trajectory &Traj,int BackRange) const{
    double rminSq=DBL_MAX;
    if(Traj.time.size()<1)
        for(int k=max(0,int(Traj.time.size()-BackRange-1));
            k<Traj.time.size() or Traj.time[k]>Traj.time[0]+tMin;
            k++)
        {
            rminSq=min(rminSq,Traj.uVec[k].normSqu());
            if(rminSq<rsqApproach)return true;
        }
    return false;
}

bool FieldTailor::stopCriterion(Trajectory & Traj, int &IStop) const {

    if(timeExceeded(Traj)){
        IStop=Traj.time.size()-1;
        return true;
    }

    if(Traj.time.back()<Traj.time[0]+tMin){
        IStop=Traj.time.size();
        return false;
    }

    // pass zero
    for(;IStop<Traj.time.size();IStop++){
        if(   not Traj.uVec[IStop  ].inGoing(Pulse::current.vecA(Traj.time[IStop ],false) )
              and Traj.uVec[IStop-1].inGoing(Pulse::current.vecA(Traj.time[IStop-1],false))){
            Traj.time.resize(IStop+1);
            Traj.uVec.resize(IStop+1);
            IStop=Traj.time.size()-1;
            return true;
        }
    }
    return false;
}

FieldTailor::FreeMotionInDipole::FreeMotionInDipole(Pulse* Pul,unsigned int Order):OdeStep("dum",0),pulse(Pul),_order(Order)
{
    _order=std::max(4,int(_order+_order%2)); // even, and at least 4
    OrthogonalLegendre().quadrature(_order/2,quad,weig);
    for(int k=0;k<quad.size();k++){
        quad[k]=(quad[k]+1.)/2.;
        weig[k]=weig[k]/2.;
    }
}

PhaseSpace & FieldTailor::FreeMotionInDipole::step(PhaseSpace &Vec, double Tstart, double Tstep)
{
    // remember, A[0]=Az,A[1]=Ax,A[2]=Ay
    for(int k=0;k<quad.size();k++){
        VectorReal a=pulse->vecA(Tstart+quad[k]*Tstep,false);
        for (int l=0;l<3;l++)Vec[l]-=a[l]*(weig[k]*Tstep);
    }
    return Vec;
}

double FieldTailor::angle(double Time, VectorReal Precol) const {
    Precol-=Pulse::current.vecA(Time);
    double ang;
    ang=acos( Precol[0]/sqrt(Precol.normSqu()));
    if(Precol[1]<0)ang=2*math::pi-ang;
    return ang;
}

double FieldTailor::operator()(const column_vector& PulseParams) const {

    // set pulse to current parameters
    for(int k=0;k<PulseParams.size();k++)
        Pulse::current.updateParameter(parRange.line(k)-1,parRange.name(k),PulseParams(k));

    // get start time (=peak field)
    double tstart=startTime();

    // align field at start time with z-axis
    Pulse::current.alignFieldZ(tstart);
    Trajectory tra=runSingle(tstart,PhaseSpace());
    if(not timeExceeded(tra))zeroIn(tra);

    // compute penalty
    double res;
    res=(tra.uVec.back().q()-positionTarget).normSqu()/(rsqApproach*4.e-2);
    VectorReal p=tra.uVec.back().p();
    if(sigmaAngle>0.){
        res+=pow((angle(tra.time.back(),p)-angleTarget)/sigmaAngle*0.1,2);
    }

    if(sigmaMomentumSqu>0.){
        VectorReal A=Pulse::current.vecA(tra.time.back());
        res+=pow(((p-A).normSqu()-momentumSquTarget)/sigmaMomentumSqu,2);
    }
    return res;
}
