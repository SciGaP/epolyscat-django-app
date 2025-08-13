// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef TIMEPROPAGATOROUTPUT_H
#define TIMEPROPAGATOROUTPUT_H

#include <vector>
#include <complex>
#include <string>
#include <iostream>
#include <fstream>
#include "mpiWrapper.h"

/// derive a class from this for your own output
class Wavefunction;
class Coefficients;
//class Operator;
class OperatorAbstract;
class DiscretizationDerived;
class Timer;
class PlotCoefficients;
class Discretization;
class OperatorMap;
class ChannelsSubregion;
class CoefficientsViewDeep;
class OperatorTree;
class CoefficientsLocal;
class DiscretizationSurface;
class TimePropagatorOutput{

    static int obj_count;
    std::vector<std::vector<std::complex<double> > > coefsT;
    std::ofstream* expecStream;
    int _expecSample; // calculation of expectation values scales poorly, keep low unless needed otherwise
    unsigned int countPrint;
    unsigned int countWrite;
    std::string dir;
    const PlotCoefficients* wfPlot;
    ChannelsSubregion* _channels;
    bool _checkGrowingNorm;
    double _tStart,_tEnd;
    double _writeInterval; ///< minimal time to elaps between writes (default: = 0)
    double _lastTimeWritten; ///< this time is actually on file
    double _printInterval;  ///< print at these intervals (default: 10 times during total)
    std::vector<OperatorAbstract*> _expecOp;
    std::vector<const OperatorAbstract*> _writeOp;
    std::vector<Coefficients> temp;
    std::vector<std::ofstream *> streamBin,streamAsc;
    Timer *_timer;
    DiscretizationDerived* _discSpec;

    void openExpec(std::string ExpecFile);
public:
    TimePropagatorOutput(
            double PrintInterval=0.,      /**< Print at these exact time-intervals (0=1/10 of total time) */
            std::string WriteDir="",      /**< Wave-function save directory (""=do not save) */
            const PlotCoefficients* WfPlot=0,   /**< Plot wave function */
            std::vector<const OperatorAbstract*> Write=std::vector<const OperatorAbstract*>(0), /**< apply these operators to wavefunction and save, filename=operator name */
            double WriteInterval=-1., /**< Minimal wave function save interval (default=PrintInterval, 0=every valid time step) */
            bool WriteAscii=false    /**< write also in ascii (at WriteInterval) */
            );
    void setInterval(double TBeg, double TEnd){_tStart=TBeg;_tEnd=TEnd;} ///< interval for progress monitoring
    virtual void write(const Coefficients *C, double Time, bool Force=false); ///< save every converged valid wave function
    void coefsFFT();
    void flush() const; ///< flush all output streams
    void close(); ///< flush all output streams

    void printStart(double accuracy, double FixStep, double CutEnergy,std::string ODEname); ///< initial info about time-propagator

    double tStart(){return _tStart;}
    double tEnd(){return _tEnd;}
    virtual void print(const Coefficients *Wf, double Time, double TimeCPU); ///< default printout

    double printInterval(){return _printInterval;}
    double writeInterval(){return _writeInterval;}
    double nextWriteTime(){return _lastTimeWritten+_writeInterval;} ///< write if time>= nextWriteTime and increment to nextWriteTime=time+writeInterval
    void lastTimeWritten(double Time){_lastTimeWritten=Time;}

    void addExpec(OperatorAbstract* Op);
    void sampleExpec(int Sample){_expecSample=Sample;}
    void checkGrowingNorm(bool Check){_checkGrowingNorm=Check;}

    Timer* timer(){return _timer;}

    //HACK - transform coefficients before Fourier transforming
    TimePropagatorOutput& withChannelsSubregion(ChannelsSubregion* Channels);

    /// apply op at each Print and write to FileName, default: Op->name()]
    TimePropagatorOutput& withApplyAndWrite(const OperatorAbstract* Op, std::string FileName="");
};


#endif // TIMEPROPAGATORTimePropagatorOutput_H
