// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "timePropagator.h"

#include "timeCritical.h"
#include "coefficientsGlobal.h"
#include "index.h"
#include "discretization.h"

//resolve forward declarations
#include "wavefunction.h"
#include "printOutput.h"
#include "plot.h"
#include "parameters.h"
#include "readInput.h"
#include "odeStep.h"
#include "derivativeFlat.h"
#include "parallel.h"
#include "mpiWrapper.h"
#include "pulse.h"
#include "log.h"

using namespace std;
using namespace tools;

double TimePropagator::_cutE,TimePropagator::_applyThreshold;

TimePropagator::TimePropagator(OdeStep<LinSpaceMap<Coefficients>, Coefficients> *Ode, TimePropagatorOutput *Out, double Precision, double FixStep)
    : h(0.01),h_max(100.),h_fix(FixStep),hSum(0.),hSquareSum(0.),desiredPrecision(Precision),
      successfulTimeSteps(0),failedTimeSteps(0),out(Out),ode(Ode),odeLocal(0)
{}

TimePropagator::TimePropagator(OdeStep<LinSpaceMap<CoefficientsLocal>, CoefficientsLocal> *Ode, TimePropagatorOutput *Out, double Precision, double FixStep)
    : h(0.01),h_max(100.),h_fix(FixStep),hSum(0.),hSquareSum(0.),desiredPrecision(Precision),
      successfulTimeSteps(0),failedTimeSteps(0),out(Out),ode(0),odeLocal(Ode)
{}


TIMER(propAll,)
TIMER(propStep,)
TIMER(propGather,)
void TimePropagator::propagate_intern(Coefficients *C, double & Time, double TEnd) {
    timeCritical::suspend();
    // set flag to indicate time-critical phase
    Coefficients errVec(C->idx());
    Coefficients vecT0(C->idx());

    CoefficientsLocal localE(C->idx());
    CoefficientsLocal localT0(C->idx());
    CoefficientsGlobal writeC(C->idx());
    CoefficientsGlobal * globalC=CoefficientsGlobal::view(C);
    CoefficientsLocal * localC  =CoefficientsLocal::view(globalC);
    timeCritical::setOn();

    //propagate from tStart to tEnd
    if(h_fix>h_max)ABORT("cannot fix time step at values > max time step");
    if(h_fix!=0)h=h_fix;
    if(odeLocal!=0)MPIwrapper::Bcast(globalC->storageData(),globalC->size(),MPIwrapper::master());
    vector<double> steps,errs;
    out->timer()->start();
    timeCritical::setOn();


    while(Time<TEnd) {
        STARTDEBUG(propAll);
        h=min(h,h_max);

        if(out->timer())out->timer()->monitor(Time, "h="+std::to_string(h));

        // we want to hit the final time exactly
        double hTemp=h;
        if(Time+h>TEnd-h*1.e-2)hTemp=TEnd-Time;
        h=min(h,h_max);
        if(h_fix!=0.){
            STARTDEBUG(propStep);
            if(ode!=0)ode->step(*C,Time,hTemp);
            else odeLocal->step(*localC,Time,hTemp);
            STOPDEBUG(propStep);
            successfulTimeSteps++;
            Time+=hTemp;
            Time=Threads::max(Time);
            STARTDEBUG(propGather);
            if(ode!=0)writeC=*C;
            else Parallel::allGather(&writeC,localC);
            STOPDEBUG(propGather);
            out->write(&writeC,Time);
            hSum+=h;
            hSquareSum+=h*h;
        }
        else {
            if(ode!=0)vecT0=*C;
            else localT0=*localC;

            STARTDEBUG(propStep);
            if(ode!=0)ode->stepError(*C,Time, hTemp,errVec);
            else odeLocal->stepError(*localC,Time,hTemp,localE);
            STOPDEBUG(propStep);

            double maxNrm;
            bool accept;
            if(ode!=0)accept=ode->acceptStep(desiredPrecision,h,h,maxNrm=Threads::max(errVec.norm()));
            else accept=odeLocal->acceptStep(desiredPrecision,h,h,maxNrm=Threads::max(localE.norm()));

            if(accept){
                successfulTimeSteps+=1;
                Time+=hTemp;
                STARTDEBUG(propGather);
                if(abs(Time-TEnd)<abs(h*1.e-10))Time=TEnd;
                if(ode!=0)writeC=*C;
                else Parallel::allGather(&writeC,localC);
                STOPDEBUG(propGather);
                out->write(&writeC,Time);
                hSum+=h;
                hSquareSum+=h*h;
            } else {
                failedTimeSteps++;
                if(ode!=0)*C=vecT0;
                else *localC=localT0;
            }


            if(h<1e-9 and TEnd-Time>1.e-8){
                // several tiny steps - terminate
                steps.push_back(h);
                errs.push_back(maxNrm);
                if(steps.size()>5){
                    string mess;
                    mess="step size underflow:\nSteps"+tools::str(steps)
                            +"\nError"+tools::str(errs)+"\nTimes "+tools::str(Time)+" "+tools::str(TEnd);
                    ABORT(mess);
                }
            }
            else {
                steps.clear();
                errs.clear();
            }
        }
        STOPDEBUG(propAll);
    }
    timeCritical::suspend();

    double killTime;
    ReadInput::main.read("DEBUG","killTime",killTime,tools::str(DBL_MAX),"kill if time larger than this",1,"DEBUGkill");
    if(killTime<Time)DEVABORT(Sstr+"killed on DEBUGkill="+killTime);

    if(ode==0)*C=writeC;
    out->timer()->stopTimer();

    // force write at tEnd
    if(out->nextWriteTime()<TEnd+(1.-1.e12)*out->writeInterval()){
        out->lastTimeWritten(TEnd-out->writeInterval());
        out->write(C,Time,true);
    }

    timeCritical::resume();
}

string TimePropagator::info() const {
    string s;
    if(ode!=0)s+=ode->info();
    else s+=odeLocal->info();
    if(h_fix==0.){
        s+=", accept/reject "+tools::str(successfulTimeSteps)+"/"+tools::str(failedTimeSteps);
        // with step size control we have double steps
        s+=", average step="+tools::str(hSum/successfulTimeSteps/2,2);
        s+=", step variance="+tools::str(sqrt(hSquareSum/successfulTimeSteps-pow(hSum/successfulTimeSteps,2))/2,2);
    }
    else
        s+=", fixed step="+tools::str(h_fix);
    return s;
}

void TimePropagator::print()
{
    if(desiredPrecision<1.e-11)PrintOutput::warning(Str("Very high accuracy of")+desiredPrecision+"requested, recommended values >= 1e-11");
    if(desiredPrecision>1.e-6) PrintOutput::warning(Str("Very low accuracy of")+desiredPrecision+"requested, recommended values <= 1e-6");

    if(_applyThreshold>0.){
        PrintOutput::lineItem("Threshold for operator application",_applyThreshold);
        PrintOutput::paragraph();
        if(desiredPrecision<_applyThreshold*0.09)
            PrintOutput::warning("High threshold for operator application - recommended: threshold < 10*accuracy");
    }
}

void TimePropagator::read(ReadInput &Inp, double &tBeg, double &tEnd, double &tPrint, double &tStore, double &accuracy, double &cutE, double & FixStep,
                          double & ApplyThreshold, std::string &Method)
{

    double tPulseEnd=-DBL_MAX,tPulseBeg=DBL_MAX;
    if(Pulse::current.apotMax()!=0.){
        tPulseBeg=Pulse::gettBegin();
        tPulseEnd=Pulse::gettEnd();
        tBeg=tPulseBeg;
        tEnd=tPulseEnd;
    }
    Inp.texdocuCategoryAdd("TimePropagation","begin,end,print,store,accuracy,cutEnergy,fixStep,operatorThreshold,method",
                           R"tex(
                           All parameters controlling time-propagation.
                           )tex","08,10,12,13");
    Inp.read("TimePropagation","begin",tBeg,tools::str(tBeg,12),"begin time of propagation")
            .texdocu(R"tex(
                     Start time of propagation. Defaults to beginning of laser pulse, if a pulse was defined.
                     )tex");
    Inp.read("TimePropagation","end",  tEnd  ,tools::str(tEnd,12),"end time of propagation")
            .texdocu(R"tex(
                     End time of propagation. Defaults to end of laser pulse, if a pulse was defined.
                     {\bf Caution:}
                     For spectra calculation without iSurff, time-propagation should be several optical cycles after the end of the  laser pulse,
                     Total duration determines achievable spectral resolution, in particular low-energy spectral accuracy.
                     Determine suitable durations by convergence studies.
                     )tex");
    ;
    if(tBeg<-DBL_MAX*1.e-2 or (tEnd>DBL_MAX*1.e-2)){
        PrintOutput::warning(Str("very large time-integration range: [")+tBeg+","+tEnd+"]");
        if(tBeg<-DBL_MAX/2 or (tEnd>DBL_MAX/2))ABORT("near-infinite time-interval. Maybe need to specify TimePropagator:begin and/or end?");
    }
    Inp.read("TimePropagation","print",tPrint,tools::str((tEnd-tBeg)/20,4),"printout intervals")
            .texdocu(R"tex(
                     Time interval between prints (default: 1/20 of total propagation time).
                     Time-propagation will be split into intervals of length \lcode{print}.
                     At the beginning of propagation and after each interval screen output and large output files,
                     such as \nameref{docu:Output:wavefunction} and checkpoints
                     will be saved. Do not choose this too small is it will cost compute time and disc space.
                     )tex");
    if(tStore<0.)tStore=tPrint;
    Inp.read("TimePropagation","store",tStore,tools::str(tStore,14),"minimal interval for file output - default=print or sample interval for spectrum,0=every valid step")
            .texdocu(R"tex(
                     Fine-grained saves without screen output. As propagation is not in general with equidistant time-steps,
                     \lcode{store} defines the minimal time interval after which the result of the next valid time-step is saved.
                     Quantities saved at this are expectation values and \lcode{tSurff} surfaces.
                     For spectral resolution $E_{\text{max}}$ you need \[\text{store}\gtrsim 2\pi/E_{\text{max}}\].
                     )tex");
    Inp.read("TimePropagation","accuracy",accuracy,"1.e-8","accuracy control")
            .texdocu(R"tex(
                     Approximate accuracy of the wave function amplitude.
                     In absence of absorption you should see the wave function norm constant on that scale.\\
                     {\bf Note:} you find info on average step size and its variance at the end of the output file.
                     If variance is not very large, consider a fixed step calculation. For choice of step size, see \nameref{docu:TimePropagation:cutEnergy}.
                     Caution: not rigourous measure, rather handle to control the accuracy. Verify needed accuracy by convergence studies. \\
                     Caution: only propagation accuracy, no meaning for basis set truncation errors.

                     )tex");
    Inp.read("TimePropagation","cutEnergy",cutE,tools::str(DBL_MAX),"remove energies above this")
            .texdocu(R"tex(
                     If specified, subspace-projections $Q$ will be applied during time-propagation.
                     \[
                     \ddt |\Psi\r = -i Q \mH Q |\Psi\r, Q=\sum_q |q\r\l q|
                     \]
                     where
                     \[
                     \mP |q\r = |q\r E_q, E_q< E_{\text{cut}}
                     \]
                     $\mP$ is an operator defined in Operator:\nameref{docu:Operator:projection}, defaults to
                     $H_0$ Operator:\nameref{docu:Operator:hamiltonian}.
                     \\Application examples: \nameref{docu:tutorial:10},\nameref{docu:tutorial:20}
                     )tex");
    Inp.read("TimePropagation","fixStep",FixStep,"0.","fixed step size (0 - automatic control)")
            .texdocu(R"tex(
                     Use a fixed step-size instead of automatic step size control. If applicable, can gain up to a factor 2 in speed.
                     Infer required step size from output, see \nameref{docu:TimePropagation:accuracy}.
                     When using a \nameref{docu:TimePropagation:cutEnergy} of $E_{\text{cut}}$, \lcode{fixStep} $\propto (E_{\text{cut}})\inv1$.
                     See application in \nameref{docu:tutorial:12},\nameref{docu:tutorial:23}.
                     )tex");
    ApplyThreshold=-1.;
    Inp.read("TimePropagation","operatorThreshold",ApplyThreshold,"-1","consider block=0, if result remains below this treshold")
            .texdocu(R"tex(
                     Will skip operation applications that contribute less than this threshold to $\ddt\Psi$.
                     Typically $\approx 10 \times$\lcode{accuracy}. Can speed up calculations, e.g. when initiall very few angular momenta are involved.
                     \\{\bf Usage:}\nameref{docu:tutorial:11}
                     \\Caution: always crosscheck against selecte calculations without that constraint.
                     \\See application in \nameref{docu:tutorial:13}.
                     )tex");
    Inp.read("TimePropagation","method",Method,"RK4","choices: RK4...classical RK, Butcher67...order 6 Butcher, CrankNicolson...(experimental)")
            .texdocu(R"tex(
                     Use of RK4 is recommended, no relevant gains were observed with Butcher67.
                     )tex");

    if(Pulse::current.apotMax()!=0.){
        double propEps=0.00001*Units::convert(1.,"OptCyc");
        if(tPulseBeg<tBeg-propEps or tEnd+propEps<tPulseEnd)
            PrintOutput::warning("time propagation starts or ends during pulse: prop=["
                                 +tools::str(tBeg)+","+tools::str(tEnd)+"], pulse=["
                                 +tools::str(tPulseBeg)+","+tools::str(tPulseEnd)+"]");
    }

    // for display - otherwise handed back as through reference
    _cutE=cutE;
    _applyThreshold=ApplyThreshold;
}

void TimePropagator::propagate(Wavefunction* Wf, double tEnd, const std::string Mode) {
    timeCritical::setOff();

    // initial wave function may have been set up only on master
    // VERY SUBTLE TRAP: if view of coefs has already been created, must not re-arrange storage
    if(Wf->coefs->storageData()==0)Wf->coefs->treeOrderStorage();
    MPIwrapper::Bcast(Wf->coefs->data(),Wf->coefs->size(),MPIwrapper::master());

    // propagation backwards in time
    if(Wf->time>tEnd){
        DEVABORT(Sstr+"backward propagation disabled, end:"+tEnd+"wf:"+Wf->time);
        //        backProp = true;
        //        Wf->time = -Wf->time;
        //        tEnd = -tEnd;
    }

    // if output has no interval, set it here
    // assign internal pointers and intermediate storage
    double nextTime=Wf->time,printStep=out->printInterval();
    if(printStep==0.)printStep=(tEnd-out->tStart())/10.;
    if(printStep<0)ABORT("cannot do backward time step: h="+str(printStep));
    if(h_fix!=0. and
            (remainder(tEnd-Wf->time,h_fix)>(tEnd-Wf->time)*1.e-10 or
             remainder(printStep,h_fix)>(tEnd-Wf->time)*1.e-10))PrintOutput::DEVwarning
            ("propagation time or print interval are not integer multiple of fixed time step - saved times will not be equidistant",1);

    if(Mode.find("Start")!=string::npos){
        out->printStart(desiredPrecision,h_fix,_cutE,ode?ode->name():odeLocal->name());
        out->print(Wf->coefs,Wf->time,0.);
    }
    timeCritical::setOn();
    while(nextTime<tEnd) {
        nextTime=min(nextTime+printStep,tEnd);
        if(abs(nextTime-tEnd)<printStep*1.e-3)nextTime=tEnd;
        propagate_intern(Wf->coefs,Wf->time,nextTime);
        LOG_MEM("propagation");

        //print some simple output to screen and file
        out->print(Wf->coefs,nextTime,out->timer()->secs());
    }
    timeCritical::setOff();

    if(Mode.find("Stop")!=string::npos){
        out->timer()->monitorDone(Wf->time);

        // dump FFT if so desired (hackish)
        out->coefsFFT();
        if(MPIwrapper::isMaster())
            out->close();
    }
    timeCritical::suspend();
}

string TimePropagator::str(unsigned int Brief) const{
    string s,b;
    s="tolerance=";
    b=tools::str(desiredPrecision,1);
    s+=b;
    if(Brief>0)return b;
    return s;
}

void TimePropagator::finalize() const{
    if(out){
        out->flush();
        out->close();
    }
}
