// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "timePropagatorOutput.h"

#include <complex>
#include "threads.h"
#include "timer.h"
#include "timeCritical.h"
#include "checkpoint.h"
#include "wavefunction.h"
#include "index.h"
//#include "operator.h"
#include "printOutput.h"
#include "plot.h"
#include "readInput.h"
#include "coefficientsWriter.h"
#include "channelsSubregion.h"

#include "operatorDefinition.h"

#include "parallelOperator.h"
#include "timePropagatorOutput.h"
#include "coefficientsViewDeep.h"
#include "discretizationtsurffspectra.h"
#include "discretizationSurface.h"
#include "timePropagator.h"
#include "eigenTools.h"

#ifdef _USE_FFTW_
#include <fftw3.h>
#endif
using namespace std;

TIMER(outWrite,)
TIMER(outWrite1,)
TIMER(outWrite2,)
TIMER(outWrite3,)
TIMER(outWrite4,)
TIMER(outWrite5,)

int TimePropagatorOutput::obj_count = 0;

bool checkForGrowingNorm=false;

void TimePropagatorOutput::openExpec(std::string ExpecFile){
    if(expecStream)return; // is already open
    if(not MPIwrapper::isMaster())return; // only master writes expec
    Checkpoint chPt(ExpecFile.substr(0,ExpecFile.rfind("/")));

    if(ExpecFile!=""){
        // only master writes expectation values - do not create empty files
        if(chPt())std::rename((ExpecFile).c_str(),(ExpecFile+"_crashed").c_str());
        expecStream=new ofstream();
        expecStream->open(ExpecFile.c_str(),(ios_base::openmode) ios::beg);
        ifstream ifil((ExpecFile+"_crashed").c_str(),(ios_base::openmode) ios::beg);
        if(ifil){
            std::string line;
            while(getline(ifil,line)){
                if(line[0]!='#'){
                    std::stringstream slin(line);
                    double time;
                    slin>>time;
                    if(time+_printInterval*1.e-10>chPt.time())break;
                }
                *expecStream<<line<<std::endl;
            }
        }
        std::remove((ExpecFile+"_crashed").c_str());
    }
}


TIMER(propIntern,)
TimePropagatorOutput::TimePropagatorOutput(
        double PrintInterval,std::string WriteDir, const PlotCoefficients *WfPlot, std::vector<const OperatorAbstract *> Write,
                                           double WriteInterval, bool WriteAscii)
    :expecStream(0),_expecSample(100),countPrint(0),countWrite(0),dir(WriteDir),wfPlot(WfPlot),_channels(0),_checkGrowingNorm(true),
      _tStart(DBL_MAX),_tEnd(-DBL_MAX),_writeInterval(WriteInterval),_printInterval(PrintInterval),_writeOp(Write),_discSpec(0)
{

    obj_count++;
    _timer=new Timer("propIntern"+tools::str(obj_count),"","timePropagator");
#ifndef _DEVELOP_
    _timer->setShow(0);
#endif
    if(wfPlot!=0){
        double t=_printInterval;
        if(Units::isDefined("OptCyc"))t=Units::convert(t,"DEFAULT_SYSTEM","OptCyc");
        wfPlot->setPlotInterval(t,false);
    }

    if(WriteDir!="" and MPIwrapper::isMaster()){
        Checkpoint chPt(WriteDir);
        for(unsigned int k=0;k<_writeOp.size();k++) {
            const Index* jI=_writeOp[k]->iIndex;
            Coefficients cK(jI);
            const Coefficients* jK=Threads::join(cK);
            if (Threads::isMaster()){
                std::string fileName = _writeOp[k]->name;
                if(chPt())std::rename((WriteDir+fileName).c_str(),(WriteDir+fileName+"_crashed").c_str());

                streamBin.push_back(new ofstream((WriteDir+fileName).c_str(),(ios_base::openmode) ios::beg|ios::binary));
                jK->idx()->write(*streamBin.back()); // coefficients outputs go with Index info

                if(chPt()){
                    // copy crashed into new up to checkpoint time
                    ifstream ifil((WriteDir+fileName+"_crashed").c_str(),(ios_base::openmode) ios::beg|ios::binary);
                    if(ifil){
                        Coefficients fC(jK->idx());
                        bool head=true;
                        double timeOnFile;
                        do{
                            if (not fC.read(ifil,head))break;
                            tools::read(ifil,timeOnFile);
                            if(timeOnFile>=chPt.time())break;
                            fC.write(*streamBin.back(),false);
                            tools::write(*streamBin.back(),timeOnFile);
                            head=false;
                        } while(true);
                        std::remove((WriteDir+fileName+"_crashed").c_str());
                    }
                    else {
                        streamBin.pop_back();
                        PrintOutput::warning(Sstr+"has checkpoint at"+chPt.time()+"but no matching output file"+(WriteDir+fileName)+" - file will not be written");
                    }
                }
                if(WriteAscii or ReadInput::main.flag("DEBUGascii","force ascii output for all data files")){
                    if(chPt.time()==-DBL_MAX)
                        streamAsc.push_back(new ofstream((WriteDir+_writeOp[k]->name+"Ascii").c_str(),(ios_base::openmode) ios::beg));
                    else
                        PrintOutput::warning("cannot use ASCII output when restarting from checkpoint",1);
                }
            }
        }
    }

    // storage for converting functions before write
    for(unsigned int k=0;k<_writeOp.size();k++){
        temp.push_back(Coefficients(_writeOp[k]->iIndex));
        temp.back().treeOrderStorage();
    }

    // always write first function
    _lastTimeWritten=-DBL_MAX;

    // minimal time-interval for writing
    if(WriteInterval<0.)_writeInterval=_printInterval;

}

TimePropagatorOutput & TimePropagatorOutput::withApplyAndWrite(const OperatorAbstract *Op, std::string FileName){
    _writeOp.push_back(Op);
    temp.push_back(Coefficients(_writeOp.back()->iIndex));

    if(MPIwrapper::isMaster()){
        if(FileName=="")FileName=Op->name;
        Coefficients cK(Op->iIndex);
        const Coefficients* jK=Threads::join(cK);
        if (Threads::isMaster()){
            streamBin.push_back(new ofstream((dir+FileName).c_str(),(ios_base::openmode) ios::beg|ios::binary));
            if(not streamBin.back()->is_open())ABORT("could not open output file "+dir+FileName);
            jK->idx()->write(*streamBin.back());
        }
    }
    return *this;
}

bool inClosedInterval(double X,double From,double To,double Eps=1.e-10){
    if(Eps==0)return From<=X and X<=To;
    return From-Eps*max(abs(From),1.)<X and X<To+Eps*max(abs(To),1.);
}
void TimePropagatorOutput::printStart(double accuracy, double FixStep, double CutEnergy, std::string ODEname)
 {
     if(_tStart==0. and _tEnd==0.)return;

     PrintOutput::title("PROPAGATION PARAMETERS");
     PrintOutput::paragraph();

     PrintOutput::newRow();
     PrintOutput::rowItem("begin");
     PrintOutput::rowItem("end");
     PrintOutput::rowItem("print");
     PrintOutput::rowItem("store");
     PrintOutput::rowItem("fixStep");
     PrintOutput::rowItem("accuracy");
     PrintOutput::rowItem("cutEnergy");
     PrintOutput::rowItem("method");

     PrintOutput::newRow();
     PrintOutput::rowItem(_tStart);
     PrintOutput::rowItem(_tEnd);
     PrintOutput::rowItem(_printInterval);
     if(_expecSample!=1)PrintOutput::rowItem(tools::str(_writeInterval,2)+"["+tools::str(_expecSample)+"]");
     else if(_writeInterval==0)PrintOutput::rowItem("all");
     else                      PrintOutput::rowItem(_writeInterval);

     if(FixStep>0.){
         PrintOutput::rowItem(FixStep);
         PrintOutput::rowItem("no control");
     } else {
         PrintOutput::rowItem("variable");
         PrintOutput::rowItem(accuracy);
     }
     PrintOutput::rowItem(CutEnergy);
     PrintOutput::rowItem(ODEname);

     PrintOutput::paragraph();

     if(_expecSample>1 and _writeInterval>0)
        PrintOutput::message(Sstr+"Expectation values sampled only every"+_expecSample+"'th store time");

}

void TimePropagatorOutput::print(const Coefficients *Wf, double Time, double TimeCPU){

    timeCritical::suspend();
    double growingNorm=1.;
    double time=Time;
    if(_printInterval!=DBL_MAX){

        if(countPrint==0){
            PrintOutput::title("TIME PROPAGATION ("+ReadInput::main.output()+")");


            if(Time>_tStart){
                PrintOutput::paragraph();
                PrintOutput::message(Sstr+"resuming from time"+Time+"after start at"+_tStart);
            }

            if(_expecOp.size()>2 or (_expecOp.size()==2 and _expecOp[1]->name!="H0")){
                // list definitions of expectation value operators
                PrintOutput::paragraph();
                for(size_t k=_expecOp.size()==2?1:2;k<_expecOp.size();k++){
                    PrintOutput::lineItem("<"+_expecOp[k]->name+">",_expecOp[k]->def()+" at hierarchy "+_expecOp[k]->idx()->hierarchy());
                    PrintOutput::newLine();
                }
            }

            if(_expecOp.size()>1){
                PrintOutput::paragraph();
                PrintOutput::lineItem("Initial <"+_expecOp[1]->name+">",_expecOp[1]->matrixElementUnscaled(*Wf,*Wf).real(),"",14);
                PrintOutput::newLine();
            }

            PrintOutput::paragraph();
            PrintOutput::newRow();
            PrintOutput::rowItem("   CPU ");
            PrintOutput::rowItem(" (%)");
            PrintOutput::rowItem("   Time");
            for(unsigned int i=0;i<_expecOp.size();i++){
                if(_expecOp[i]->name.length()>16)
                    PrintOutput::rowItem("<"+_expecOp[i]->name.substr(0,13)+"...>");
                else
                    PrintOutput::rowItem("<"+_expecOp[i]->name+">");
            }
            //HACK re-establish distribution pattern of operators
            for(size_t k=0;k<_expecOp.size();k++){
                ParallelOperator::setDistribution(_expecOp[k]);
            }
            for(auto o: _expecOp){
                o->update(Time,Wf);
                if(o->def().find("MeanEE")!=string::npos)dynamic_cast<OperatorTree*>(o)->updateNonLin(Time,Wf);
            }
        }
        PrintOutput::newRow();
        PrintOutput::rowItem(TimeCPU);
        PrintOutput::rowItem((Time-_tStart)/(_tEnd-_tStart)*100.);

        // write near-zero time as 0
        if(abs(time)<abs(_tStart)*1.e-12)time=0.;
        PrintOutput::rowItem(time,6);


        double complexEps=1.e-6;
#ifdef _DEVELOP_
        complexEps=1.e-10; // more sensitive for developer
#endif
        vector<unsigned int>compExp;
        for(unsigned int i=0;i<_expecOp.size();i++){
            _expecOp[i]->update(Time,Wf);
            complex<double> expec=_expecOp[i]->matrixElementUnscaled(*Wf,*Wf);
            expec=Threads::sum(expec);
            PrintOutput::rowItem(real(expec),14);
            if(abs(expec.imag())>abs(expec.real())*complexEps  and abs(expec.imag())>complexEps)compExp.push_back(i);
            if(_expecOp[i]->name.find("Overlap")!=string::npos)growingNorm=abs(expec);
        }
        for(unsigned int k=0;k<compExp.size();k++)
            if(_expecOp[compExp[k]]->name.substr(0,3)!="Ovr"){
                PrintOutput::warning(Sstr+"complex expectation value of operator "
                                     +_expecOp[compExp[k]]->name+"(threshold  ="+complexEps+")",5);
            }

        double t=time;
        if(Units::isDefined("OptCyc"))t=Units::convert(time,"DEFAULT_SYSTEM","OptCyc");
        vector<string> head(1,string("time="+tools::str(time,8)));
        const Coefficients *joinedWf=Threads::join(*const_cast<Coefficients*>(Wf));
        // only plot if either standard parallel Wf or joinedWf on master of Threads
        if(wfPlot and  (joinedWf==Wf or Threads::isMaster()))
            wfPlot->plot(*joinedWf,dir+wfPlot->briefName()+tools::str(t,5),head);

        PrintOutput::flush();
        if(growingNorm-1.>1.e-7)checkForGrowingNorm=true;
        if(Threads::isMaster())Checkpoint::write(dir,Time,joinedWf);

        write(Wf,Time,true); // force write whenever there is print
        flush();
        //        if(growingNorm>1.e5)MPIwrapper::Abort(1);

        countPrint++;
        timeCritical::resume();
    }
}

void TimePropagatorOutput::flush() const {
    for(auto b: streamBin)b->flush();
    for(auto a: streamAsc)a->flush();
    if(expecStream)expecStream->flush();
    PrintOutput::flush();
}

void TimePropagatorOutput::close(){
    if(expecStream)expecStream->flush();
    for(unsigned int c=0;c<streamBin.size();c++){
        streamBin[c]->close();
        if(streamBin[c]) delete streamBin[c];
        streamBin[c] = 0;
    }
    for(unsigned int c=0;c<streamAsc.size();c++){
        streamAsc[c]->close();
        if(streamAsc[c]) delete streamAsc[c];
        streamAsc[c] = 0;
    }
}

void TimePropagatorOutput::write(const Coefficients *C, double Time, bool Force){

    if(_lastTimeWritten>=Time)return; // never write backward or duplicate times
    if(_channels)_channels->average(C, Time);
    if(not Force and Time<nextWriteTime())return;
    STARTDEBUG(outWrite);
    _lastTimeWritten=Time;

    // header on expec file
    if(expecStream and not expecStream->tellp()){
        for(unsigned int k=0;k<_expecOp.size();k++)*expecStream<<"# "+_expecOp[k]->name<<endl;
        *expecStream<<"#      Time            CPU ";
        for(unsigned int i=0;i<_expecOp.size();i++)*expecStream<<setw(20)<<"<"+_expecOp[i]->name+">";
        *expecStream<<endl;
    }

    double t=Units::isDefined("OptCyc")? Units::convert(Time,"DEFAULT_SYSTEM","OptCyc") : Time;
    if(wfPlot){
        // on this level, true plots will only be done if in append mode
        const Plot* plot=dynamic_cast<const Plot*>(wfPlot);
        if(not plot or plot->isAppend()){
            const Coefficients *joinedC=Threads::join(*const_cast<Coefficients*>(C));
            // only Threads master writes
            if(Threads::isMaster())
                wfPlot->plot(*joinedC,dir+wfPlot->briefName(),std::vector<std::string>(),tools::str(t,12));
        }
    }

    // apply operators and write wave functions
    STARTDEBUG(outWrite2);
    for(unsigned int k=0;k<_writeOp.size();k++){
        const_cast<OperatorAbstract*>(_writeOp[k])->axpy(1.,*C,0.,temp[k],Time);
        const Coefficients *joinedTemp=Threads::join(temp[k]);
        if(Threads::isMaster()){
            joinedTemp->write(*streamBin[k],false);
            tools::write(*streamBin[k],Time);
            if(streamAsc.size()>k){
                *streamAsc[k]<<Time;
                joinedTemp->print(*streamAsc[k]);
                streamAsc[k]->flush();
            }
        }
    }
    STOPDEBUG(outWrite2);

    // expectation values (if specified)
    STARTDEBUG(outWrite1);
    if((Force or (_expecSample and not (countWrite%_expecSample))) and _expecOp.size()){
        if(expecStream)*expecStream<<setw(20)<<setprecision(12)<<Time<<setw(12)<<setprecision(4)<<_timer->secs();
        for(unsigned int i=0;i<_expecOp.size();i++){
            STARTDEBUG(outWrite3);
            complex<double> expec=_expecOp[i]->matrixElementUnscaled(*C,*C);
            STOPDEBUG(outWrite3);
            STARTDEBUG(outWrite4);
            MPIwrapper::Barrier();
            if(i>1 and abs(expec.imag())>abs(expec.real())*1.e-12  and abs(expec.imag())>1.e-12
                    and _expecOp[i]->name!="Hamiltonian")
                PrintOutput::warning("complex expectation value of operator "+_expecOp[i]->name,5);
            double expecRe=Threads::sum(expec.real());
            STOPDEBUG(outWrite4);
            STARTDEBUG(outWrite5);
            if(expecStream)*expecStream<<setw(20)<<setprecision(12)<<expecRe;
            STOPDEBUG(outWrite5);
        }
        if(expecStream)*expecStream<<std::endl;
    }
    if(_checkGrowingNorm and C->idx()->overlap()){
        double growingNorm=C->idx()->overlap()->matrixElementUnscaled(*C,*C).real();
        growingNorm=Threads::sum(growingNorm);
        //        MPIwrapper::Bcast(&growingNorm,1,MPIwrapper::master());
        if(growingNorm>1.2){
            PrintOutput::warning(Str("norm>1, possible instability - try smaller scaling angle or smaller operator threshold"),1,0,
                                 "Possible reasons - trouble  shoot: \
                                 \n - scaling angle too large - use smaller\
                                 \n - TimePropagator:accurcy or fixStep too large \
                                 \n - basis in absorption range too small - use larger \
                                 \n - polynomial order too high - use lower \
                                 \n - spectral cut in time-propagator fails - increase cut energy \
                                 "
                                 );
        }
        if(growingNorm<=1.)checkForGrowingNorm=false;
    }
    STOPDEBUG(outWrite1);
    countWrite++;
    STOPDEBUG(outWrite);
}


void TimePropagatorOutput::addExpec(OperatorAbstract* Op){
    _expecOp.push_back(Op);
    openExpec(ReadInput::main.output()+"expec");
}



void TimePropagatorOutput::coefsFFT() {
    if(not OperatorAbstract::flat)return;
    if(coefsT.size()==0)return;
#ifdef _USE_FFTW_
    ofstream fftStream((ReadInput::main.output()+"coefs").c_str(),(ios_base::openmode) ios::beg);
    PrintOutput::paragraph();
    unsigned int nOmega=coefsT[0].size()/400;
    PrintOutput::message("Fourier transform of Coefficients on file "+ReadInput::main.output()+"coefs");

    fftw_complex *in, *out;
    fftw_plan p=0;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * coefsT[0].size());
    out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * coefsT[0].size());
    vector<double> smooth(coefsT[0].size());
    for  (unsigned int k=0;k<smooth.size();k++)
        smooth[k]=pow(sin(k*math::pi/double(smooth.size())),8);
    for(unsigned int k=0;k<coefsT.size();k++){
        p = fftw_plan_dft_1d(coefsT[0].size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        for (unsigned int l=0;l<coefsT[k].size();l++){*(in+l)[0]=norm(coefsT[k][l])*smooth[l];*(in+l)[1]=0;}

        fftw_execute(p); /* repeat as needed */

        // output can be read into gnuplot
        double sum=0.;
        for (unsigned int l=0;l<coefsT[k].size();l++){
            sum+=(pow(*(in+l)[0],2)+pow(*(in+l)[1],2));
            if((l+1)%nOmega==0){
                fftStream<<sum<<endl;
                sum=0.;
            }
        }
        fftStream<<endl;

    }
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
#endif
}

TimePropagatorOutput& TimePropagatorOutput::withChannelsSubregion(ChannelsSubregion* Channels){
    _channels=Channels;
    return *this;
}
