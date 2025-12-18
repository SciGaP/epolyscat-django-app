// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef TIMER_H
#define TIMER_H
#include <vector>
#include <string>
#include <climits>
#include "time.h"
#include "stdio.h"
#include "abort.h"
#include "platformSpecific.h"

#ifdef _TIMER_DEBUG_
#undef _TIMER_OFF_
#endif

#define TIMERSAMPLEFACTOR 1024

// set up a file static object TIMERname and define macros for start/stop
#ifndef _TIMER_OFF_

#define TIMER(name,group)          static Timer          TIMER ## name ( #name, #group, __FILE__);
#define TIMERSAMPLE(name,group)    static TimerSample    TIMER ## name ( #name, #group, __FILE__);
#define TIMERRECURSIVE(name,group) static TimerRecursive TIMER ## name ( #name, #group, __FILE__);
#define START(name)           TIMER ## name.start();
#define STOP(name)            TIMER ## name.stopTimer();

#define TIMERBEGIN(name)      TimerList::setList( #name, __FILE__);
#define TIMEREND(name)        TimerList::unsetList( #name, __FILE__);
#define TIMERSTART(name,group) TimerList::start( #name, #group, __FILE__);
#define TIMERSTOP(name,group)  TimerList::stopTimer(#name, __FILE__);

#if defined(_TIMER_DEBUG_) || defined(_DEVELOP_)
#define STARTDEBUG(name) TIMER ## name.start();
#define STOPDEBUG(name)  TIMER ## name.stopTimer();
#else
#define STARTDEBUG(name) ;
#define STOPDEBUG(name) ;
#endif

#else

#define TIMER(name,group) ;
#define TIMERSAMPLE(name,group)   ;
#define TIMERRECURSIVE(name,group) ;
#define TIMERBEGIN(name) ;
#define TIMEREND(name) ;
#define TIMERSTART(name,group) ;
#define TIMERSTOP(name,group) ;
#define STARTDEBUG(name)       ;
#define STOPDEBUG(name)        ;
#define START(name)       ;
#define STOP(name)        ;
#define TIMERSAMPLE(name,group)    ;
#define TIMERRECURSIVE(name,group) ;
#endif

/// \ingroup Abort
class Timer
        /** @brief In-line timer using pairs of macros START(timername) ... STOP(timername) */
        /// Usage:
        ///
        ///
        /// #include "timer.h"
        ///
        /// TIMER(timername,groupname); // = static Timer TIMERtimername("timername","groupname",__FILE__);
        /// void myFunction(){
        ///   ...code...
        /// START(timername) // = TIMERtimername.start();
        ///   ...code...
        /// STOP(timername) // = TIMERtimername.stop();
        ///   ...code...

        /// The timername string needs to be different for all Timers within the file
        /// <br>Group: connect to a timer of that name in the timers defined so far, or create a new group
        /// <br>Group is optional, if omitted, no group will be assigned
        /// <br>compiler flag -D_TIMER_OFF_   deactivate all timers (re-compile ALL timed code)
        /// <br>compiler flag -D_TIMER_DEBUG_ activate STARTDEBUG/STOPDEBUG (re-compile ALL timed code)

        /// <br><br>Alternate usage (slightly more overhead):
        ///
        ///
        /// void myFunction(){
        ///   ...code...
        /// // open a timer region
        /// TIMERBEGIN(region); // = TimerList::setList("region",__FILE__);
        ///   ...code...
        /// TIMERSTART(name,group) // = TimerList::start("name","group",__FILE__);
        ///   ...code...
        /// TIMERSTOP(name) // = TimerList::stopTimer("name",__FILE__);
        ///   ...code...
        /// // close the region
        /// TIMEREND(region); // = TimerList::unsetList("region",__FILE__);

        /// There can only be one region opened at a time
        /// <br>While a region is open, pairs of TIMERSTART/TIMERSTOP can be called anywhere in the code, also in subroutines
        /// <br>If no region is open, TIMERSTART/TIMERSTOP has no effect (but causes function call overhead)
        /// <br>!!! NOTE !!!
        /// <br>Conflicts "timer logics out of synch" may arise if the sequence of timer calls in a region changes at repeated runs through the region
        /// <br>(for speed, timers are identified in sequence of their calls in the region)

{

    friend class TimerList;
    static std::vector<Timer*> table; // any timer is pointed to from this table

    // need timerOverflow%TIMERSAMPLEFACTOR=0 for sampling to work reliably in case of overflow
    static const unsigned long timerOverflow=INT_MAX/2-(INT_MAX/2)%TIMERSAMPLEFACTOR;

    inline void addDiff(timespec & sum,timespec & tstart, timespec & tend){
        if(tstart.tv_nsec<tend.tv_nsec){
            sum.tv_nsec+=tend.tv_nsec-tstart.tv_nsec;
            sum.tv_sec +=tend.tv_sec -tstart.tv_sec;
        } else {
            sum.tv_sec +=tend.tv_sec-tstart.tv_sec-1;
            sum.tv_nsec+=1000000000+tend.tv_nsec -tstart.tv_nsec;
        }
        if(sum.tv_nsec>1000000000){
            sum.tv_nsec-=1000000000;
            sum.tv_sec++;
        }
    }
    void insertTable();

protected:
    bool running;       ///< indicates whether timer is running
    timespec cpu;       ///< cpu time
    timespec wall;      ///< wall clock time
    timespec tstart;    ///< start time
    timespec tstop;     ///< stop time
    timespec tmon;      ///< most recent monitor time
    unsigned long calls; ///< number of calls to timer
    unsigned int overf; ///< overflow count
    std::string file;   ///< file from which called
    std::string name;   ///< string to associate with the timer
    std::string info;   ///< further info on the time (used by monitor)
    Timer * group;      ///< timers can come in (assumed nested) hierarchies
    std::vector<Timer*> member; ///< pointers to lower levels of hierarchy
    unsigned int recursiveLevel; ///< count level of recursive call
    int sample;  ///< timer sampled a the given period
    std::ofstream * monOut; ///< output stream for monitor output
    int _show;  ///< 1...show in writes, 0...do not show in write

    static void overhead(std::ostream &Out=std::cout); // estimate the timer overhead
    static bool compare_file_group_name(const Timer* A, const Timer * B);
    static bool compare_time(const Timer* A, const Timer * B);
    static unsigned long overLimit;
    inline void incrCall() {
        calls++;
        // fast test for imminent overflow:
        if(calls>timerOverflow){
            calls=1;
            overf+=1;
        }
    }

public:
    static Timer generalMonitor;
    /// active by default, deactivate by compiling all code with -D_TIMER_OFF_
    Timer(const std::string & Name, const std::string & Group, const std::string & File);
    static std::vector<std::string> stopAll();

    void setShow(int Show){_show=Show;} ///<set show level, 10...always show, 0...do not show
    /// start timer, set up upon first start
    virtual inline void start()
    {
        if(calls==0)insertTable();
        if(running)ABORT("cannot start timer "+file+":"+name+", is already running"+
                               "\nfor recursive use: startRecursive/stopRecursive");
        incrCall();
        running=true;
        platformSpecific::current_utc_time(&tstart);
    }

    /// stop timer, add elapsed time to total time
    virtual inline void stopTimer(){
        platformSpecific::current_utc_time(&tstop);
        addDiff(cpu,tstart,tstop);
        if(not running)ABORT("cannot stop timer "+file+":"+name+", is not running"+
                                   "\nfor recursive use: startRecursive/stopRecursive");
        running=false;
    }

    static const std::string currentDateTime();
    static void write(std::string File); ///< formatted timer summary to file
    static void write(std::ostream& Out);///< formatted timer summary to any ostream
    std::string shortWall();             ///< return 3 digit time string

    void monitor(double Parameter, const std::string Info="", const std::string Out="",bool Done=false);
    void monitorDone(double Parameter){monitor(Parameter,"","",true);}
    double secs(){return cpu.tv_sec+cpu.tv_nsec*1.e-9;}       ///< total CPU time in seconds
    double wallSecs(){return wall.tv_sec+wall.tv_nsec*1.e-9;} ///< total wall clock time in seconds
    double dSecs() const; ///< most recent difference of tstart/tstop timer pair

};

/// sample timer, for very frequent calls to very fast functions
class TimerSample : public Timer {
public:
    TimerSample(const std::string & Name,const std::string & Group, const std::string & File)
        : Timer(Name,Group,File){sample=TIMERSAMPLEFACTOR;}
    /// for very frequent calls, just sample
    inline void start(){     if(calls%TIMERSAMPLEFACTOR==0)Timer::start();else incrCall();}
    inline void stopTimer() {if(calls%TIMERSAMPLEFACTOR==1)Timer::stopTimer();}
};

/// recursive timer: called only at top-level of recursive functions
class TimerRecursive : public Timer {
public:
    TimerRecursive(const std::string & Name,const std::string & Group, const std::string & File)
        : Timer(Name,Group,File){sample=1;}
    inline void start(){recursiveLevel++;if(recursiveLevel==1)Timer::start();}
    inline void stopTimer() {if(recursiveLevel==1)Timer::stopTimer();recursiveLevel--;}

    // when this is called, it should show only one call count (although it  recuresively descends 10 levels)
    static void test(int level=0);

};

// for dynamic, yet fast timer use
class TimerList {

    static std::map<std::string,TimerList*> _timerLists;
    static TimerList* _currentList;
    static Timer * currentList(std::string Name,std::string Group);

    std::vector<Timer*> list;
    std::string name,file;
    unsigned int n;
public:
    TimerList(std::string Name, std::string File);
    static void setList(std::string ListName, std::string ListFile);
    static void unsetList(std::string ListName, std::string ListFile);
    static void start(std::string Name, std::string Group,std::string File);
    static void stopTimer(std::string Name, std::string File);
};


#endif // TIMER_H
