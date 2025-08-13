#include "farm.h"
#include "tools.h"
#include "readInput.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <set>

#include "multiIndex.h"
#include "printOutput.h"
#include "inputGenerator.h"

Farm::Farm()
    :_comFarm(MPIwrapper::communicator()),_comFlock(MPIwrapper::communicator()),_inp(""){
    for(int k=0;k<MPIwrapper::Size(_comFarm);k++)_flock.push_back(k);
    _statusTimer=new Timer("Farm","",__FILE__);
}

Farm::Farm(std::shared_ptr<InputGenerator> Generator):Farm(){
    _inputGenerator=Generator;
}

static std::set<int> running;
static double _previousStatus=-DBL_MAX/2;

bool Farm::runCompleted() const{
    return _inputGenerator->runCompleted(tools::stringStrip(_inp,ReadInput::inputCopy));
}

bool Farm::next(){
    MPIwrapper::setCommunicator(_comFarm);
    if(MPIwrapper::Size()==1  or dynamic_cast<InputSingle*>(_inputGenerator.get())!=0){
        // no slaves or trivial single input
        _inp=_inputGenerator->nextInput();
        assigned.push_back(_inp);
        if(assigned.back()=="")assigned.pop_back();
        return _inp!="";
    }

    if(MPIwrapper::isMaster(_comFarm)){

        while(""!=(_inp=_inputGenerator->nextInput())){
            // wait for report from slave
            std::string report;
            int slav=MPIwrapper::Recv(report,MPIwrapper::anySource());
            PrintOutput::message(Sstr+"assigning "+tools::stringStrip(_inp,"/"+ReadInput::inputCopy)+"to"+slav);
            // assign new task;
            assigned.push_back(_inp);
            MPIwrapper::Send(_inp,slav);
            running.insert(slav);
            if((_statusTimer->secs()-_previousStatus)>1.){
                if(_previousStatus>-DBL_MAX/2)status();
                _previousStatus=_statusTimer->secs();
            }

            std::string s;
            for(auto i: running)s+=" "+tools::str(i);
        }
        std::string done("");
        for(auto s: running){
            MPIwrapper::Send(done,s);
        }
        return false;
    }
    else {
        if(MPIwrapper::isMaster(_comFlock)){
            // flock masters reports back
            MPIwrapper::Send(_inp,0);
            // receives new assignement
            MPIwrapper::Recv(_inp,0);
        }
        // broadcast to whole flock
        MPIwrapper::setCommunicator(_comFlock);
        MPIwrapper::Bcast(_inp,MPIwrapper::master());
    }
    return _inp!="";
}

std::string Farm::inp(){
    return _inp;
}

void Farm::fork(){

    // do not fork if single task (or scalar run)
    if(MPIwrapper::Size()==1 or _inputGenerator->nInps()==1)return;
    if(not MPIwrapper::isMaster())PrintOutput::off("Screen");

    int siz=0;
    if(MPIwrapper::isMaster())siz=_inputGenerator->nInps();
    MPIwrapper::Bcast(&siz,1,MPIwrapper::master());
    int groupSize=std::max(1,std::min((MPIwrapper::Size(_comFarm)-1)/siz,_maxFlockSize));

    std::vector<std::vector<int>>_flocks(1);
    for(int k=1;k<MPIwrapper::Size(_comFarm);k++){
        _flocks.back().push_back(k);
        if(_flocks.back().size()>=(size_t)groupSize){
            // create flock on present process
            if(std::find(_flocks.back().begin(),_flocks.back().end(),MPIwrapper::Rank(_comFarm))!=_flocks.back().end()){
                _flock=_flocks.back();
            }
            if(k+groupSize<MPIwrapper::Size(_comFarm))_flocks.push_back({});
        }
    }
        _comFlock=MPIwrapper::setCommunicator(_flock,_comFarm);
}

std::string Farm::status(){
    if(dynamic_cast<InputSingle*>(_inputGenerator.get()))return "";

    if(MPIwrapper::isMaster(_comFarm)){
        std::vector<std::string> res;
        std::string head("no header found");
        for(auto inp: assigned){
            std::string s,h;
            std::string mon=tools::stringStrip(inp,ReadInput::inputCopy)+"mon";
            std::ifstream stream(mon.c_str(),std::ios::in);
            if(stream.is_open()){
                s="... no monitor file";
                std::getline(stream,h);
                std::getline(stream,s);
                if(h!="")head="   "+tools::cropString(h);
            }
            stream.close();
            if(s==""){
                s="no Status info yet";
            } else {
                s=s.substr(s.find("Stat:")+5);
                s=tools::cropString(s);
            }
            res.push_back(s);
        }

        PrintOutput::end();
        PrintOutput::title("Status of assigned runs");
        PrintOutput::newRow();
        PrintOutput::rowItem("run");
        if(res.size())head.resize(res[0].size(),' ');
        head.back()='.';
        PrintOutput::rowItem(head);
        for(size_t k=0;k<res.size();k++){
            PrintOutput::newRow();
            PrintOutput::rowItem(tools::stringStrip(assigned[k],ReadInput::inputCopy));
            PrintOutput::rowItem(res[k]);
        }
        PrintOutput::paragraph();
    }
    _inputGenerator->print();
    return "";
}
