// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "timer.h"
#include <stdio.h>
#include <algorithm>    // std::sort

#include "platformSpecific.h"

using namespace std;

#ifndef _NO_TOOLS_
#include "tools.h"
#include "printOutput.h"
using namespace tools;
#else
//=== auxiliary methods --- replacement for functions in tools.h ===============
//* a range of auxiliary methods that would be part of other includes
* here we try to keep the number of includes low and define the few function
* locally (file-static) instead
*
* static keyword:
    *    in a file (outside any class) it declares a symbol with
  *    "file scope" - it is visible ONLY to code in the same file
  *    this keeps any other code with accidentally the same name from linking to this symbol
  */
  #include <cmath>    // for pow(x,n)
  #include <iostream> // basic i/o methods
  #include <iomanip>  // control of i/o: setw(..),precision(..)
  #include <fstream>  // file i/o stream
  #include <sstream>  // string i/o stream
  static std::string fileBase(const std::string & File){
    size_t i0=File.rfind("/")+1;
    if(i0==std::string::npos)i0=0;
    return File.substr(i0,File.rfind(".")-i0);
}
static std::string lcropString(std::string s){while(s.length()!=0){if(s[           0]!=' '&&s[           0]!='\r'&&s[           0]!='\n')break;s.erase(0,           1);};return s;}
static std::string rcropString(std::string s){while(s.length()!=0){if(s[s.length()-1]!=' '&&s[s.length()-1]!='\r'&&s[s.length()-1]!='\n')break;s.erase(s.length()-1,1);};return s;}
static std::string cropString(std::string s){return lcropString(rcropString(s));}

static std::string str(int P){ostringstream oss; oss<< P;return oss.str();}

static void printWarning(std::string Message){cout<<"\n !!! "+Message+" !!!\n";}
//=== end of auxiliary methods ================================================
#endif

Timer Timer::generalMonitor("Monitor","",__FILE__);

std::vector<Timer*> Timer::table;
unsigned long Timer::overLimit=std::pow(2,30);

Timer::Timer(const string &Name, const string &Group, const string &File)
    :running(false),calls(0),overf(0),name(Name),group(0),recursiveLevel(0),monOut(0),_show(10)
{
    sample=1;
    file=tools::fileBase(File);
    cpu.tv_sec=0;
    cpu.tv_nsec=0;
    wall.tv_sec=0;
    wall.tv_nsec=0;
    tmon.tv_sec=0;
    tmon.tv_nsec=0;
    group=0;
    return;
}

TimerList* TimerList::_currentList=0;
map<std::string,TimerList*> TimerList::_timerLists;

TimerList::TimerList(string Name, string File)
    :name(Name),n(0){
    file=tools::fileBase(File);
}

/// open a list, or create new if needed (call through macro TIMERBEGIN(listname)
void TimerList::setList(std::string ListName /** unique within file */, std::string ListFile /** use __FILE__ */){
    if(_timerLists.count(ListName+ListFile)==1){
        if(_currentList!=0)ABORT("cannot setList, timer list still set: "+_currentList->name+" from file "+_currentList->file);
        _currentList=_timerLists[ListName+ListFile];
    }
    else {
        _currentList=new TimerList(ListName,ListFile);
        _timerLists[ListName+ListFile]=_currentList;
    }
    _currentList->n=0;
}

/// close list, must be called before any list can be restarted; call through macro TIMEREND(listname)
void TimerList::unsetList(string ListName, string ListFile){
    if(_currentList!=_timerLists[ListName+ListFile])ABORT("timer list not running: "+ListName+" in file "+ListFile);
    _currentList=0;
}

/// start timer in list, create new if needed; macro: TIMERSTART(name, group)
void TimerList::start(string Name,string Group, string File){
    if(_currentList==0)return;
    if(_currentList->list.size()<_currentList->n+1){
        if(_currentList->list.size()!=_currentList->n)ABORT("timer logics out of synch");
        // find timer with same name
        unsigned int nt;
        for(nt=0;nt<Timer::table.size();nt++)if(Timer::table[nt]->name==_currentList->name+"."+tools::cropString(Name) and
                                                Timer::table[nt]->file==_currentList->file+"@"+tools::fileBase(File))break;
        if(nt<Timer::table.size())_currentList->list.push_back(Timer::table[nt]);
        else                      _currentList->list.push_back(new Timer(_currentList->name+"."+tools::cropString(Name),Group,
                                                                         _currentList->file+"@"+tools::fileBase(File)));
        _currentList->list[_currentList->n]->start();
    }
    else
        _currentList->list[_currentList->n]->start();
}

/// stop timer in list; macro: TIMERSTOP(name, group)
void TimerList::stopTimer(string Name, string File){
    if(_currentList==0)return;
    _currentList->list[_currentList->n]->stopTimer();
    _currentList->n++;
}

double Timer::dSecs() const {
    double delta;
    if(tstart.tv_nsec<tstop.tv_nsec){
        delta=(tstop.tv_nsec-tstart.tv_nsec)*1.e-9+tstop.tv_sec-tstart.tv_sec;
    } else {
        delta=tstop.tv_sec-tstart.tv_sec-1+(1000000000+tstop.tv_nsec -tstart.tv_nsec)*1.e-9;
    }
    return delta;
}

const std::string Timer::currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%m/%d %X", &tstruct);
    return buf;
}

void Timer::monitor(double Parameter, const string Info, const string Out, bool Done){
    // try set up
    if(Info!="")info=Info; // update info
    if(monOut==0){
        if(Out!=""){
            monOut=new ofstream(Out.c_str(),(ios_base::openmode) ios_base::trunc);
            if(not monOut->is_open()){
                PrintOutput::warning("could not open monitor output file \""+Out+"\" --- no monitor output will be printed",10);
                delete monOut;
                monOut=0;
            }
        } else {
            PrintOutput::DEVwarning("before monitoring timer "+name+" set output file by monitor(paramter,info,monitorFile) ",1);
            monOut=0;
        }
        if(not running or monOut==0)return; // not set up or not started yet
    }

    // get current time
    timespec tend;
    platformSpecific::current_utc_time(&tend);
    if(tend.tv_sec-tmon.tv_sec<2 and not Done) return; // under minimal interval and not final call
    tmon=tend;

    if(running){
        // advance time to present (like stop/start pair)
        addDiff(cpu,tstart,tend);
        tstart=tend;
    } else {
        if(calls==0){
            PrintOutput::DEVwarning("for monitoring, must start timer "+name+" by start()");
        }
    }

    monOut->seekp(ios_base::beg); // rewind file
    *monOut<<"        Date   Time     (CPU)      Param.   \tInfo...."<<endl;
    if(Done)*monOut<<"Stat: *"<<currentDateTime();
    else    *monOut<<"Stat:  "<<currentDateTime();
    *monOut<<" ("<<setw(7)<<secs()<<")";
    *monOut<<" "<<setw(8)<<Parameter<<"\t"<<info<<endl<<flush;

    if(Done){
        monOut->close();
        delete monOut;
        monOut=0;
    }

}

void Timer::write(std::string File){
#ifndef _TIMER_OFF_
    cout<<"\n *** timing information on file "+File+" ***"<<endl;
    std::ofstream out;
    out.open(File.c_str());
    write(out);
    out.close();
#else
    cout<<"\n *** timers disabled by _TIMER_OFF_ ***"<<endl;
#endif
}
void Timer::insertTable(){
    // check for uniqueness and insert into table
    for(unsigned int n=0;n<table.size();n++){
        if(table[n]->file==file){
            if(table[n]->name==name){
                if(table[n]->group==group){
                    ABORT("duplicate timer "+file+"::"+name);
                }
            }
        }
    }
    table.push_back(this);
}

std::vector<std::string> Timer::stopAll(){
    std::vector<std::string> names;
    for(auto t: table){
        if(t->running){
            names.push_back(t->name);
            t->stopTimer();
        }
    }
    return names;
}

void Timer::write(std::ostream & Out){
#ifndef _TIMER_OFF_
    std::sort(table.begin(),table.end(),Timer::compare_time);
    Out<<"\n                     File: Name                          calls   cpu(total)   cpu/call [samp]"<<endl;
    for(unsigned int n=0;n<table.size();n++){
        Timer * t=table[n];
        if (not t->_show)continue;
        Out.width(25);
        string f=t->file;
        if(f.find("@")!=0)f=f.substr(f.find("@")+1);
        Out<<f<<": ";
        Out<<setw(25)<<left<<t->name<<right;
        if(t->calls==0)Out<<"..."<<endl;
        else{
            Out<<setw(10)<<t->calls
              << setprecision(3)
                 //              <<setw(12)<<double(t->wall.tv_sec+t->wall.tv_nsec*1.e-9)*t->sample
              <<setw(12)<<double(t->cpu.tv_sec+t->cpu.tv_nsec*1.e-9)*t->sample
             <<setw(12)<<(double(t->cpu.tv_sec+t->cpu.tv_nsec*1.e-9)*t->sample)/t->calls;
            if(t->sample>1)Out<<"  ["<<setw(4)<<t->sample<<"]";
            if(t->overf!=0)Out<<"  Calls overflow "<<t->overf<<" x "<<timerOverflow;
            Out<<endl;
        }
    }
    overhead(Out);
#else
    Out<<"\n  --- timers disabled (for enabling, remove _TIMER_OFF_ from compilation) ---\n"<<endl;
#endif
}

// sorting of timers on printout
bool Timer::compare_file_group_name(const Timer* A,const Timer * B) {
    if(A->file==B->file){ // first files first
        if(A->group==B->group)return A->name<B->name;
        if(A->group==0)return true; // non-group first
        if(B->group==0)return false;
        return A->group->name<B->group->name;
    }
    return A->file<B->file;
}
// sorting of timers on printout
bool Timer::compare_time(const Timer* A,const Timer * B) {
    if(A->cpu.tv_sec==B->cpu.tv_sec)return A->cpu.tv_nsec>B->cpu.tv_nsec;
    return A->cpu.tv_sec>B->cpu.tv_sec;
}

string Timer::shortWall(){
    if(wallSecs()>=60*60*24*7*52)return ">1y";
    int secs=wallSecs();
    if(secs<0)ABORT("cannot convert negative time: "+str(secs));
    if(secs<10       )return " "+str(secs)+"s";
    if(secs<60       )return str(secs)+"s";
    if(secs<60*10    )return str(secs/60)+"m"+str(secs%60/6);
    if(secs<60*100   )return str(secs/60)+"m";
    if(secs<60*60*10 )return str(secs/(60*60))+"h"+str(secs%(60*60)/(60*6));
    if(secs<60*60*100)return str(secs/(60*60))+"h";
    if(secs<60*60*24*10)return  str(secs/(60*60*24))+"d"+str(int(secs%(60*60*24)/(60*60*2.4)));
    if(secs<60*60*24*100)return  str(secs/(60*60*24))+"d";
    if(secs<60*60*24*7*10)return  str(secs/(60*60*24*7))+"w"+str(int(secs%(60*60*24*7)/(60*60*24*0.7)));
    if(secs<60*60*24*7*52)return  str(secs/(60*60*24*7))+"w";
    return ">1y";
}

void Timer::overhead(std::ostream &Out){
    Timer overFull("overFull","",__FILE__);
    Timer overSamp("overSample","",__FILE__);
    Timer dum("dum","",__FILE__);
    TimerSample sam("sam","",__FILE__);
    const int N=100000;

    overFull.start();
    for(int n=0;n<N;n++){
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
        dum.start();
        dum.stopTimer();
    }
    overFull.stopTimer();
    Out<<"timer overhead/call (sec): "<<overFull.secs()/(16*N)<<'\n';

    overSamp.start();
    for(int n=0;n<N;n++){
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
        sam.start();
        sam.stopTimer();
    }
    overSamp.stopTimer();
    Out<<"timer overhead/call (sec): "<<overSamp.secs()/(16*N)<<" (sample="<<sam.sample<<")\n";
    // remove temporary timers from table
    Timer::table.pop_back();
    Timer::table.pop_back();
    Timer::table.pop_back();
    Timer::table.pop_back();
}

TIMERRECURSIVE(recursive,)
void TimerRecursive::test(int level){
    START(recursive);
    if(level==0)cout<<"testRecursive timer: ";
    cout<<level<<" ";
    if(level>10)return;
    level++;
    TimerRecursive::test(level);
    STOP(recursive);
    if(level==0)cout<<endl;
}

