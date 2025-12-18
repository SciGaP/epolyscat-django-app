// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "threads.h"

#include "mpiWrapper.h"
#include "readInput.h"
#include "index.h"
#include "coefficients.h"
#include "parallel.h"
#include "basisGridQuad.h"
#include "basisSub.h"
#include "timeCritical.h"
#include "readInput.h"
#include "debugInfo.h"


#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <map>


static std::map<std::string, const Index*> _parent;
static std::map<std::string, const Index*> _jointIndices;
static std::map<std::string, Threads*> _threads;
static std::set<std::string> _isSet;
static std::map<std::string, std::string> _threadsLevel;

static MPI_Comm allThreads=0;
static MPI_Comm singleThread=0;
static MPI_Comm thisThread=0;

void Threads::setup(MPI_Comm Communicator){
    if(allThreads)return; // already set up
    singleThread=MPIwrapper::setCommunicator({MPIwrapper::Rank(Communicator)},Communicator);
    thisThread=allThreads=Communicator;
    MPIwrapper::setCommunicator(thisThread);
}
MPI_Comm Threads::all(){return allThreads;}
MPI_Comm Threads::single(){return singleThread;}

double Threads::max(double Val){
    if(MPIwrapper::communicator()!=Threads::single())return Val;
    MPIwrapper::setCommunicator(all());
    MPIwrapper::AllreduceMAX(&Val,1);
    MPIwrapper::setCommunicator(single());
    return Val;
}

std::complex<double> Threads::sum(std::complex<double> Val){
    if(MPIwrapper::communicator()!=Threads::single())return Val;
    MPIwrapper::setCommunicator(all());
    MPIwrapper::AllreduceSUM(&Val,1);
    MPIwrapper::setCommunicator(single());
    return Val;
}
double Threads::sum(double Val){
    if(MPIwrapper::communicator()!=Threads::single())return Val;
    MPIwrapper::setCommunicator(all());
    MPIwrapper::AllreduceSUM(&Val,1);
    MPIwrapper::setCommunicator(single());
    return Val;
}

int Threads::rank() {
    return MPIwrapper::Rank(all());
}

bool Threads::isMaster(){
    return MPIwrapper::isMaster(Threads::all());
}

int determineRows(const Index* Tdx,const Index* Idx){
    const Index* fdx=Idx->firstFloor();
    const Index* tdx=fdx->axisIndex(Tdx->axisName().substr(7));
    if(not tdx or tdx->depth()>fdx->depth()+1)DEVABORT(Tdx->str()+"\njoin only at floor level or one below");
//    DEVABORT(Tdx->str()+"\njoin only at floor level or one below");
    return fdx==tdx ? 0 : fdx->basis()->size();
}

Threads::Threads(const Index *Tdx){
    _T.reset(new Coefficients(Tdx));
    _C.reset(new Coefficients(join(Tdx)));
    _rows=determineRows(_T->idx(),_C->idx());
}

bool Threads::set(const Coefficients* C){
    if(_isSet.count(C->idx()->hash()))return true;
    if(not isThread(C->idx()))return false;

    timeCritical::suspend();
    _isSet.insert(C->idx()->hash());
    Threads* thr=new Threads(C);
    if(Threads::isMaster()){
        _threads[C->idx()->hash()]=thr;
    }
    else
        delete thr;
    timeCritical::resume();

    return true;
}

Index* assemble(const Index* Idx){

    if(MPIwrapper::communicator()!=Threads::single())return 0;
    MPIwrapper::setCommunicator(Threads::all());

    Index *tdx=0;
    if(not MPIwrapper::isMaster()){
        // write into string
        std::ostringstream idxFile;
        Idx->write(idxFile);
        idxFile.flush(); //do we need this?
        // send string to master
        int siz=idxFile.str().length();
        MPIwrapper::Send(siz,MPIwrapper::master(),0);
        MPIwrapper::Send(const_cast<char*>(&idxFile.str()[0]),idxFile.str().length(),0,1);
    }

    if(MPIwrapper::isMaster()){
        tdx = new Index();
        tdx->childAdd(new Index(*Idx));       // first branch: sub-index on master
        for(int k=1;k<MPIwrapper::Size();k++) // remaining branches
        {   // receive string
            int siz;
            MPIwrapper::Recv(siz,k,0);
            std::string buf(siz,' ');
            MPIwrapper::Recv(&buf[0],siz,k,1);
            // read-construct from string
            std::istringstream idxFile(buf);
            tdx->childAdd(new Index(idxFile));
        }
        tdx->sizeCompute();
    }

    MPIwrapper::setCommunicator(Threads::single());
    return tdx;
}

std::string threadsLevel(const Index* Idx){
    if(MPIwrapper::communicator()!=Threads::single())return "";
    if(not _threadsLevel.count(Idx->hash()))
    {
        Index * tdx=assemble(Idx);
        if(tdx){
            if(tdx->childSize()==1){
                // (temporary for debug) single thread with trivial subbasis
                for(Index *cdx=tdx->child(0);cdx!=0;cdx=cdx->descend()){
                    if(dynamic_cast<const BasisSub*>(cdx->basis()) and
                            cdx->basis()->size()==BasisSub::superBas(cdx->basis())->size()){
                        _threadsLevel[Idx->hash()]=cdx->axisName();
                        break;
                    }
                }
            }
            else {
                // find level where bases disagree
                for(Index *ix0=tdx->child(0),*ix1=tdx->child(1);ix0!=0;ix0=ix0->descend(),ix1=ix1->descend()){
                    if(not (ix0->basis()==ix1->basis())
                            and not (*ix0->basis()==*ix1->basis())){
                        _threadsLevel[Idx->hash()]=ix0->axisName();
                        break;
                    }
                }
            }
        }
        else {
            _threadsLevel[Idx->hash()]="";
        }
        delete tdx;
        MPIwrapper::setCommunicator(Threads::all());
        MPIwrapper::Bcast(_threadsLevel[Idx->hash()],MPIwrapper::master());
        MPIwrapper::setCommunicator(Threads::single());
    }
    return _threadsLevel[Idx->hash()];
}

bool Threads::isThread(const Index *Idx){
    return threadsLevel(Idx)!="";
}

Threads::Threads(const Coefficients *C){

    std::string lev=threadsLevel(C->idx());
    if(lev=="")DEVABORT(C->idx()->str()+"\nError: no Threads level identified");


    Index * tdx=assemble(C->idx());
    if(tdx){
        tdx->setAxisName("Threads"+lev);
        if(tdx->axisName()=="Threads")ABORT(C->idx()->str()+"\nError: no Threads level identified");
        _threadedI.reset(tdx);
        _joinedI.reset(join(_threadedI.get()));
        _T.reset(new Coefficients(_threadedI.get()));
        _C.reset(new Coefficients(_joinedI.get()));
        _rows=determineRows(_T->idx(),_C->idx());
    }
}

void Threads::join(Coefficients* J, std::vector<const Coefficients*> & VF,int Rows){
    if(J->isLeaf()){
        if(Rows==0){
            // vertically merge row-wise stored matrices (easy)
            int kJ=0;
            for(const Coefficients* f: VF){
                for(size_t kf=0;kf<f->size();kf++,kJ++)
                    J->floorData()[kJ]=f->floorData()[kf];
            }
        }
        else {
            // horizontally concatenate row-wise stored matrices
            for (int row=0,kJ=0;row<Rows;row++){
                for(const Coefficients* f: VF)
                    for(size_t  kf=0;kf<f->size()/Rows;kf++,kJ++){
                        J->floorData()[kJ]=f->floorData()[row*f->size()/Rows+kf];
                    }
            }
        }
    }

    for(size_t  k=0;k<J->childSize();k++){
        std::vector<const Coefficients*> vF;
        for(const Coefficients* f: VF)vF.push_back(f->child(k));
        join(J->child(k),vF,Rows);
    }
}

void Threads::scatter(const Coefficients *J, std::vector<Coefficients *> &VF, int Rows){
    if(J->isLeaf()){
        if(Rows==0){
            // vertically split row-wise stored matrices (easy)
            int kJ=0;
            for(Coefficients* f: VF){
                for(size_t  kf=0;kf<f->size();kf++,kJ++)
                    f->floorData()[kf]=J->floorData()[kJ];
            }
        }
        else {
            // horizontally split row-wise stored matrices
            for (int row=0,kJ=0;row<Rows;row++){
                for(Coefficients* f: VF)
                    for(size_t  kf=0;kf<f->size()/Rows;kf++,kJ++){
                        //                        J->floorData()[kJ]=f->floorData()[row*f->size()/Rows+kf];
                        f->floorData()[row*f->size()/Rows+kf]=J->floorData()[kJ];
                    }
            }
        }
    }

    for(size_t  k=0;k<J->childSize();k++){
        std::vector<Coefficients*> vF;
        for(Coefficients* f: VF)vF.push_back(f->child(k));
        scatter(J->child(k),vF,Rows);
    }
}

Coefficients * Threads::join(Coefficients &C){
    if(MPIwrapper::communicator()!=Threads::single())return const_cast<Coefficients*>(&C);

    bool boolSet=Threads::set(&C);
    if(not boolSet)return const_cast<Coefficients*>(&C);; // create threads (if needed)

    // join threads
    Coefficients* jC=0;
    std::vector<const Coefficients*> vF;
    MPIwrapper::setCommunicator(Threads::all());
    if(not Threads::isMaster()){
        MPIwrapper::Send(C.data(),C.size(),MPIwrapper::master());
    } else {
        Threads* thr=_threads[C.idx()->hash()];
        *thr->_T->child(0)=C;
        vF.push_back(thr->_T->child(0));
        for (int k=1;k<MPIwrapper::Size();k++){
            MPIwrapper::Recv(thr->_T->child(k)->orderedData(),thr->_T->child(k)->size(),k);
            vF.push_back(thr->_T->child(k));
        }
        join(thr->_C.get(),vF,thr->_rows);
        jC=thr->_C.get();
    }
    MPIwrapper::setCommunicator(Threads::single());
    return jC;
}

void Threads::scatter(const Coefficients *C, Coefficients &Scattered){
    if(not Threads::set(&Scattered)){
        Scattered=*C;
        return;
    }

    // distribute to threads
    MPIwrapper::setCommunicator(Threads::all());
    if(Threads::isMaster()){
        Threads * thr=_threads[Scattered.idx()->hash()];
        std::vector<Coefficients*>vT;
        for(size_t  k=0;k<thr->_T->childSize();k++)vT.push_back(thr->_T->child(k));
        // scatter joint into threaded data
        scatter(C,vT,thr->_rows);
        Scattered=*vT[0];
        for (int k=1;k<MPIwrapper::Size();k++){
            MPIwrapper::Send(vT[k]->orderedData(),vT[k]->size(),k);
        }
    }
    else{
        MPIwrapper::Recv(Scattered.data(),Scattered.size(),MPIwrapper::master());
    }
    MPIwrapper::setCommunicator(Threads::single());
}

const Index* Threads::join(const Index *Idx){
    if(not _jointIndices.count(_parent[Idx->hash()]->hash())){
        if(Idx->axisName().substr(0,7)!="Threads")DEVABORT(Sstr+"not a Thread index and not in _jointIndices\n"+Idx->str());
        Index* jdx=new Index(*Idx->child(0));
        for(Index* tdx=jdx->axisIndex(Idx->axisName().substr(7));tdx;tdx=tdx->nodeRight()){
            tdx->setBasis(BasisSub::superBas(tdx->basis()));
            if(not tdx->isBottom())
                for(size_t  k=tdx->childSize();k<tdx->basis()->size();k++)
                    tdx->childAdd(new Index(*tdx->child(0)));
        }
        jdx->sizeCompute();
        _jointIndices[Idx->hash()]=jdx;
    }
    return _jointIndices[Idx->hash()];
}

Index *Threads::detach(const Index *Idx){
    // no detaching on single thread (may be admitted for debugging)
    if(MPIwrapper::Size(MPIwrapper::communicator())==1)return new Index(*Idx);

    if(Idx->axisName().substr(0,7)!="Threads")DEVABORT("cannot detach from non-Thread node: "+Idx->hierarchy());
    Index* detIdx=new Index(*Idx->child(MPIwrapper::Rank()));
    _parent[detIdx->hash()]=Idx;
    MPIwrapper::setCommunicator(single());
    return detIdx;
}


Index *Threads::fork(const Index *Idx, std::string AxisName){
    // no forking on single thread (may be admitted for debugging)
    if(MPIwrapper::Size(MPIwrapper::communicator())==1)return new Index(*Idx);

    Index* topIdx=new Index();
    topIdx->setBasis(BasisAbstract::factory("Vector:"+tools::str(MPIwrapper::Size())));
    topIdx->setAxisName("Threads"+AxisName);
    for(size_t n=0;n<topIdx->basis()->size();n++){
        Index* idx=new Index(*Idx);
        idx=idx->axisIndex(AxisName);
        const BasisAbstract* fBas=idx->basis();
        int k0=n*fBas->size()/MPIwrapper::Size();
        int k1=std::min(size_t((n+1)*fBas->size()/MPIwrapper::Size()),size_t(fBas->size()));
        std::vector<int> subset;
        for(int k=k0;k<k1;k++)subset.push_back(k);
        for(Index* jdx=idx->axisIndex(AxisName);jdx;jdx=jdx->nodeRight()){
            if(not (*jdx->basis()==*fBas))DEVABORT("for forking, need all bases on level equal");
            jdx->setBasis(BasisAbstract::factory(BasisSub::strDefinition(jdx->basis(),subset)));
            if(not jdx->isLeaf())
                for(size_t k=0;k<fBas->size()-subset.size();k++)jdx->childPop();
        }

        // index ownerships
        idx=idx->root();
        Parallel::setIndexOwner(idx,n);
        for(Index * jdx=const_cast<Index*>(idx->firstFloor());jdx;jdx=jdx->nodeRight())Parallel::setIndexOwner(jdx,n);

        idx->resetFloor(Idx->firstFloor()->depth());
        idx->sizeCompute();
        topIdx->childAdd(idx); // attach to main
    }
    _jointIndices[topIdx->hash()]=Idx;
    topIdx->sizeCompute();
    topIdx->resetFloor(Idx->firstFloor()->depth()+1);
    _threads[topIdx->hash()]=new Threads(topIdx);
    return topIdx;
}
