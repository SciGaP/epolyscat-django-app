// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "parallelProcess.h"

#include <map>

#include "readInput.h"
#include "parallelGrain.h"
#include "parallel.h"
#include "index.h"
#include "operatorFloor.h"
#include "coefficientsFloor.h"
#include "mpiWrapper.h"
#include "timer.h"

#include "debugInfo.h"

using namespace std;

bool ParallelProcess::debugTimer=false;

ParallelProcess::~ParallelProcess(){
    //    for(auto p: temp)delete p;
    //    for(auto p: localB)delete p;
    //    for(auto p: sendB)delete p;
    //    for(auto p: recvB)delete p;
}

ParallelProcess::ParallelProcess(unsigned int Nproc, unsigned int Numb):numb(Numb),_used(false)
{
    sendB.resize(Nproc);
    recvB.resize(Nproc);
    sendC.resize(Nproc);
    recvC.resize(Nproc);
}

double ParallelProcess::load() const {
    double l=0;
    for(unsigned int k=0;k<temp.size();k++)l+=temp[k]->load();
    for(unsigned int k=0;k<localB.size();k++)l+=localB[k]->load();
    for(unsigned int n=0;n<sendB.size();n++)for(unsigned int k=0;k<sendB[n].size();k++)l+=sendB[n][k]->load();
    for(unsigned int n=0;n<recvB.size();n++)for(unsigned int k=0;k<recvB[n].size();k++)l+=recvB[n][k]->load();
    return l;
}

void ParallelProcess::addCross(ParallelCross *Cross){
    // bypass check if thread
    if(Parallel::owner(Cross->index())!=number()){
        DEVABORT(Sstr+MPIwrapper::Rank()+Cross->index()->hash()+"Cross does not belong to process"+number());
    }
    for(unsigned int n=0;n<Cross->colBlock.size();n++)temp.push_back(Cross->colBlock[n]);
    for(unsigned int n=0;n<Cross->rowBlock.size();n++)temp.push_back(Cross->rowBlock[n]);
}

bool ParallelProcess::unused() const {return not _used;}

void ParallelProcess::setSendRecv(const string & Assign){
    //NOTE: for now names of "Assign" are mis-used
    _used=temp.size()>0;
    for(unsigned int n=0;n<temp.size();n++)
    {
        const DerivativeBlock* b=temp[n];
        if(std::find(_outFloors.begin(),_outFloors.end(),b->cInOut[1])==_outFloors.end())
            _outFloors.push_back(b->cInOut[1]);

        // get owners
        unsigned int iOwner=Parallel::owner(b->cInOut[1]->idx());
        unsigned int jOwner=Parallel::owner(b->cInOut[0]->idx());
        // define as local, if not assigned yet
        // NOTE: it is unclear why ownership would be set to Parallel::all (maybe a leftover from earlier implementation of projectors)
        if(iOwner==Parallel::none){
            if(Assign=="either")   Parallel::setIndexOwner(b->cInOut[1]->idx(),iOwner=jOwner);
            else if(Assign=="send")Parallel::setIndexOwner(b->cInOut[1]->idx(),iOwner=Parallel::all);
        }
        if(jOwner==Parallel::none){
            if(Assign=="either")      Parallel::setIndexOwner(b->cInOut[0]->idx(),jOwner=iOwner);
            else if(Assign=="receive")Parallel::setIndexOwner(b->cInOut[1]->idx(),jOwner=Parallel::all);
        }
        if(     (iOwner==number() and (jOwner==number() or jOwner==Parallel::all or jOwner==Parallel::thread)) or
                (jOwner==number() and (iOwner==number() or iOwner==Parallel::all or iOwner==Parallel::thread))){
            localB.push_back(b);
        }
        else if (iOwner==number())
            recvB[jOwner].push_back(b); // recv rhs before application
        else if (jOwner==number())
            sendB[iOwner].push_back(b); // send lhs after application
        else {
            ABORT("neither index is owned by process: "
                  +tools::str(number())+": "+tools::str(iOwner)+" <-- "+tools::str(jOwner)
                  +" - programming error: ");
        }
    }
    temp.clear();
}

TIMER(debugISend,)
void ParallelProcess::sendTo(unsigned int Recipient, std::vector<Coefficients *> &C, MPIwrapper::Buffer &Buf){
    if(Buf.size()==0)return; // nothing to be sent
    // copy into send buffer
    unsigned int pos=0;
    for(unsigned int k=0;k<C.size();k++)
        for(unsigned int l=0;l<C[k]->size();l++,pos++)Buf.val[pos]=*(C[k]->floorData()+l);
    if(pos!=Buf.val.size())MPIwrapper::Out()<<"buffer sizes do not match "<<pos<<" vs, "<<Buf.val.size()<<"\n";
    if(debugTimer)STARTDEBUG(debugISend);
    MPIwrapper::ISend(Buf.val.data(),Buf.size(),Recipient,Buf.req);
    if(debugTimer)STOPDEBUG(debugISend);
    if(Buf.req==0)ABORT(tools::str(MPIwrapper::Rank())+"request not set properly");

}

TIMER(debugRecv,)
void ParallelProcess::recvFrom(unsigned int Sender, std::vector<Coefficients *> &C, MPIwrapper::Buffer &Buf){
    if(Buf.size()==0)return;
    // receive
    if(debugTimer)STARTDEBUG(debugRecv);
    MPIwrapper::Recv(Buf.val.data(),Buf.size(),Sender);
    if(debugTimer)STOPDEBUG(debugRecv);
    // copy from receive buffer
    unsigned int pos=0;
    for(unsigned int k=0;k<C.size();k++)
        for(unsigned int l=0;l<C[k]->size();l++,pos++)*(C[k]->floorData()+l)=Buf.val[pos];
    if(pos!=Buf.val.size())MPIwrapper::Out()<<"buffer sizes do not match "<<pos<<" vs, "<<Buf.val.size()<<"\n";
}

void ParallelProcess::addFrom(unsigned int Sender, std::vector<Coefficients *> &C, MPIwrapper::Buffer &Buf){
    if(Buf.size()==0)return;
    // confirm receive
    if(debugTimer)STARTDEBUG(debugRecv);
    MPIwrapper::Recv(Buf.val.data(),Buf.size(),Sender);
    if(debugTimer)STOPDEBUG(debugRecv);
    // copy from receive buffer
    unsigned int pos=0;
    for(unsigned int k=0;k<C.size();k++)
        for(unsigned int l=0;l<C[k]->size();l++,pos++)*(C[k]->floorData()+l)+=Buf.val[pos];
    if(pos!=Buf.val.size())MPIwrapper::Out()<<"buffer sizes do not match "<<pos<<" vs, "<<Buf.val.size()<<"\n";
}

void ParallelProcess::setBuffers(const Parallel * Par){
    recvBuf.resize(recvB.size());
    for(unsigned int n=0;n<recvB.size();n++){
        unsigned int siz=0;
        for(unsigned int k=0;k<recvB[n].size();k++){
            if(std::find(recvC[n].begin(),recvC[n].end(),recvB[n][k]->cInOut[0])==recvC[n].end()){
                recvC[n].push_back(recvB[n][k]->cInOut[0]);
                siz+=recvC[n].back()->size();
            }
        }
        recvBuf[n].val.resize(siz);
    }

    sendBuf.resize(sendB.size());
    for(unsigned int n=0;n<sendB.size();n++){
        unsigned int siz=0;
        for(unsigned int k=0;k<sendB[n].size();k++){
            if(std::find(sendC[n].begin(),sendC[n].end(),sendB[n][k]->cInOut[1])==sendC[n].end()){
                sendC[n].push_back(sendB[n][k]->cInOut[1]);
                siz+=sendC[n].back()->size();
            }
        }
        sendBuf[n].val.resize(siz);
    }
}

string ParallelProcess::str() const {

    // info on traffic from/to process
    // (from)->[this]->(to), recv/send=#data{#floors}/#data{#floors}

    string s;

    unsigned int sum=0;
    sum=localB.size();
    for(unsigned int m=0;m<sendB.size();m++)sum+=sendB[m].size();
    for(unsigned int m=0;m<recvB.size();m++)sum+=recvB[m].size();
    if(sum>0)s+="blocks= "+tools::str(sum);

    string from,to;
    unsigned int toSize=0,fromSize=0,fromF=0,toF=0;
    for(unsigned int n=0;n<sendB.size();n++){
        if(sendB[n].size()>0)to+=tools::str(n)+",";
        if(sendC[n].size()>0)toF+=sendC[n].size();
        toSize+=sendBuf[n].size();
    }

    for(unsigned int n=0;n<recvB.size();n++){
        if(recvB[n].size()>0)from+=tools::str(n)+",";
        if(recvC[n].size()>0)fromF+=recvC[n].size();
        fromSize+=recvBuf[n].size();
    }

    // remove trailing comma, add buffer sizes
    string Com;
    if(from.length()>0){
        from.resize(from.length()-1);
        Com="\tRecv= "+tools::str(fromSize,0)+"{"+tools::str(fromF,0)+"}";
    }
    if(to.length()>0){
        to.resize(to.length()-1);
        if(Com!=""){
            Com.insert(5,string("/Send"));
            Com+="/"+tools::str(toSize,0)+"{"+tools::str(toF,0)+"}";
        }
        else Com="\tSend= "+tools::str(toSize,0)+"{"+tools::str(toF,0)+"}";
    }

    s+=", ("+from+")->["+tools::str(number())+"]->("+to+")"+Com;

    if(localB.size()==0 and from.length()==0 and to.length()==0)return "(unused)";

    return s;

}

void ParallelProcess::check(){
    for(unsigned int k=0;k<localB.size();k++)
        if(     Parallel::_indexOwner[localB[k]->cInOut[1]->idx()->hash()]!=number() or
                Parallel::_indexOwner[localB[k]->cInOut[0]->idx()->hash()]!=number())
            ABORT("block is not diagonal");

    for(unsigned int n=0;n<sendB.size();n++)
        for(unsigned int k=0;k<sendB[n].size();k++){
            if(     Parallel::_indexOwner[sendB[n][k]->cInOut[1]->idx()->hash()]!=n or
                    Parallel::_indexOwner[sendB[n][k]->cInOut[0]->idx()->hash()]!=number())
                ABORT("host or recipient do not match");
        }

    for(unsigned int n=0;n<recvB.size();n++){
        for(unsigned int k=0;k<recvB[n].size();k++){
            if(     Parallel::_indexOwner[recvB[n][k]->cInOut[0]->idx()->hash()]!=n or
                    Parallel::_indexOwner[recvB[n][k]->cInOut[1]->idx()->hash()]!=number())
                ABORT("host or sender do not match");
        }
    }
}
