// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "parallel.h"

#define _PARALLEL_

#include "readInput.h"
#include "discretization.h"
#include "index.h"
#include "coefficients.h"
#include "coefficientsGlobal.h"
////#include "operator.h"
#include "operatorFloor.h"
#include "derivativeBlock.h"
#include "parallelProcess.h"
#include "parallelCross.h"
#include "parallelLayout.h"
#include "operatorTensor.h"
#include "basisMat1D.h"
#include "basisMat2D.h"
#include "basisMatMulti.h"

#include "printOutput.h"

#include "mpiWrapper.h"
#include "parallelBuffer.h"
#include "parallelOperator.h"
#include "parallelLayout.h"

#include "log.h"

#include "toolsPrint.h"
#include "debugInfo.h"

bool Parallel::debugTimer=false;

using namespace std;


unsigned int Parallel::sizeNode=4;
unsigned int Parallel::sizeBoard=2;
unsigned int Parallel::sizeCPU=2;

map<const Index*,vector<unsigned int> > Parallel::indexSort;
map<std::string,unsigned int> Parallel::_indexOwner;
static std::map<std::string,const Index*> _ownedIndex;

map<std::string,unsigned int> Parallel::_operatorHost;
map<std::string,unsigned int> Parallel::_grainHost; ///< index of process that hosts iIndex,jIndex block
std::unique_ptr<ParallelLayout> Parallel::_layout;

Parallel::~Parallel(){
    //    delete proc;
    //    proc=0;
    for(auto p: cross)delete p;
    for(auto p: grain)delete p;
}

void Parallel::setSort(const Discretization *Disc){
    if(Disc==0)return;
    if(Disc->idx()->sizeCompute()==0)return; // empty discretization, no sorting

    // create given sorting
    vector<unsigned int> s;
    for(unsigned int k=0;k<Disc->idx()->firstFloor()->depth();k++)s.push_back(k);

    //HACK - will fail for hybrid hierarchies
    int itop=0;
    for(const Index* idx=Disc->idx();idx!=0;idx=idx->descend()){
        if(idx->continuity()!=Index::npos)std::swap(s[itop++],s[idx->depth()]);
    }

    indexSort[Disc->idx()]=s;

    _layout.reset(new ParallelLayout(Disc->idx(),""));
    _layout->setFloorHosts(Disc->idx());
}

void Parallel::construct(unsigned int NProc, unsigned int Level)
{
    if(Level==0){
        for(unsigned int n=0;n<NProc;n++)childAdd(new Parallel(new ParallelProcess(MPIwrapper::Size(),n)));
        return;
    }

    unsigned int nsub=0;
    if(Level==3)nsub=sizeNode;
    if(Level==2)nsub=sizeBoard;
    if(Level==1)nsub=sizeCPU;

    while(NProc>nsub){
        childAdd(new Parallel(nsub,Level-1));
        NProc-=nsub;
    }
    childAdd(new Parallel(NProc,Level-1));

    if(Level==3){
        for(Parallel* p=descend(4);p!=0;p=p->nodeRight()){
            _process.push_back(p->proc);
            _process.back()->numb=_process.size()-1;
        }
    }
}

void Parallel::setToZeroLHS() {
    for(Coefficients* f: _process[MPIwrapper::Rank()]->_outFloors)f->setToZero();
}


TIMER(debugWait,)
TIMER(debugApply,)
void Parallel::apply(const std::complex<double> &Alfa) {

    if(debugTimer)STARTDEBUG(debugApply);
    proc=process(MPIwrapper::Rank());

    unsigned int curBegin=0,curEnd=MPIwrapper::Size();
    curBegin=MPIwrapper::Rank();
    curEnd=curBegin+1;

    ParallelProcess * curProc;

    // initiate send to other processes
    for(unsigned int cur=curBegin;cur<curEnd;cur++){
        curProc=process(cur);
        for(int n=0;n<MPIwrapper::Size();n++){
            if(process(n)->recvBuf[curProc->number()].size()>0){
                if(curProc==proc)
                    proc->sendTo(process(n)->number(),process(n)->recvC[curProc->number()],process(n)->recvBuf[curProc->number()]);
            }
        }
    }

    // apply to local data and initiate send
    for(unsigned int cur=curBegin;cur<curEnd;cur++){
        curProc=process(cur);
        for(unsigned int n=0;n<curProc->sendB.size();n++){
            if(curProc->sendBuf.size()>0){
                // apply and send
                applySubset(curProc->sendB[n],Alfa);
                if(curProc==proc)
                    proc->sendTo(process(n)->number(),curProc->sendC[process(n)->number()],curProc->sendBuf[process(n)->number()]);
            }
        }
    }

    // all strictly local operations
    for(unsigned int cur=curBegin;cur<curEnd;cur++){
        curProc=process(cur);
        applySubset(curProc->localB,Alfa);
    }


    // receive rhs and apply, result is local
    for(unsigned int cur=curBegin;cur<curEnd;cur++){
        curProc=process(cur);
        for(int n=0;n<MPIwrapper::Size();n++){
            if(curProc->recvB[n].size()>0){
                // receive and apply
                if(curProc==proc)
                    curProc->recvFrom(process(n)->number(),curProc->recvC[process(n)->number()],curProc->recvBuf[process(n)->number()]);
                applySubset(curProc->recvB[n],Alfa);
            }
        }
    }

    // receive results from others
    for(unsigned int cur=curBegin;cur<curEnd;cur++){
        curProc=process(cur);
        for(unsigned int n=0;n<curProc->sendB.size();n++){
            if(curProc==proc)
                curProc->addFrom(process(n)->number(),process(n)->sendC[proc->number()],process(n)->sendBuf[proc->number()]);
        }
    }

    // collect all local request and wait for completion
    vector<MPI_Request>req;
    for(int k=0;k<MPIwrapper::Size();k++){
        if(proc->sendBuf[k].size()>0                   )req.push_back(proc->sendBuf[k].req);
        if(process(k)->recvBuf[proc->number()].size()>0)req.push_back(process(k)->recvBuf[proc->number()].req);
    }
    if(debugTimer)STARTDEBUG(debugWait);
    MPIwrapper::Waitall(req);
    if(debugTimer)STOPDEBUG(debugWait);
    if(debugTimer)STOPDEBUG(debugApply);
}

void Parallel::keepLocal(Coefficients* C){
    for(Coefficients* f=const_cast<Coefficients*>(C->firstLeaf());f!=0;f=f->nextLeaf())
        if(int(owner(f->idx()))!=MPIwrapper::Rank())f->setToZero();
}


#ifdef _OPENMP
void OMP_storage(const std::vector<const DerivativeBlock *> &Block){
    if(Block.back()->_alt==0){
        Block.back()->_alt=Block.back()->cInOut[1]->floorData();
        for(auto b=Block.begin();b<Block.end()-1;b++)
            if(Block.back()->_alt==(*b)->cInOut[1]->floorData())
                Block.back()->_alt=new std::complex<double>[Block.back()->cInOut[1]->size()];
    }
    if((Block.back())->_alt!=(Block.back())->cInOut[1]->floorData())
        memset(Block.back()->_alt,0,Block.back()->cInOut[1]->size()*sizeof(*(Block.back())->_alt));
}
#endif

TIMER(applySubset,)
void Parallel::applySubset(const std::vector<const DerivativeBlock *> &Block,std::complex<double>Alfa) const
{
    STARTDEBUG(applySubset)
            std::vector<const DerivativeBlock*>::const_iterator b;
    complex<double> alfa;
#pragma omp parallel default(shared) private(b,alfa)
    {
#pragma omp for schedule(dynamic,1) nowait
        for(b=Block.begin();b<Block.end();b++)
        {
            alfa=Alfa;
            if((*b)->oLeaf->floor()->factor()!=0){
                if(*(*b)->oLeaf->floor()->factor()==0.)continue;
            }

            // application threshold defined, skip if below
#ifndef _OPENMP
            if((*b)->eps>=0 and max(abs(alfa.real()),abs(alfa.imag()))*(*b)->cInOut[0]->norm() <= (*b)->eps)continue;
#else
            if((*b)->_skip=((*b)->eps>=0 and max(abs(alfa.real()),abs(alfa.imag()))*(*b)->cInOut[0]->norm()<=(*b)->eps))continue;
            if((*b)->_alt==0)OMP_storage(std::vector<const DerivativeBlock *>(Block.begin(),b+1));
            memset((*b)->_alt,0,(*b)->cInOut[1]->size()*sizeof(*(Block.back())->_alt));
#endif

            //NOTE: use of complex multiplications seems to be faster then rearranging real and imag parts
#ifndef _OPENMP
            (*b)->oLeaf->floor()->apply(alfa,(*b)->cInOut[0]->floorData(),(*b)->cInOut[0]->size(),
                    1.                      ,(*b)->cInOut[1]->floorData(),(*b)->cInOut[1]->size());
#else
            (*b)->oLeaf->floor()->apply(alfa,(*b)->cInOut[0]->floorData(),(*b)->cInOut[0]->size(),
                    1.                      ,(*b)->_alt,                  (*b)->cInOut[1]->size());
#endif
        }
    }
#ifdef _OPENMP
    // add from temporary storage into vector
    for(b=Block.begin();b<Block.end();b++){
        if((*b)->_skip)continue;
        if((*b)->_alt!=(*b)->cInOut[1]->floorData()){
            for(std::complex<double>*c=(*b)->cInOut[1]->floorData(),*a=(*b)->_alt;
                a!=(*b)->_alt+(*b)->cInOut[1]->size();a++,c++)*c+=*a;
        }
    }
#endif
    STOPDEBUG(applySubset)
}

string Parallel::str() const {
    string s="total load ="+tools::str(load())+"     (pre-receive)->[oper]->(post-send)    dataSize{#blocks}";
    for(int k=0;k<MPIwrapper::Size();k++){
        const ParallelProcess * procHere=process(k);
        s+="\nProcess no "+tools::str(procHere->number());
        if(procHere==0)
            s+=" (unused)";
        else
            s+=",\tloadDelta="+tools::str(int(100*(MPIwrapper::Size()*procHere->load()/load()-1.)),3)
                    +"%,\t"+procHere->str();
    }
    return s;
}

string Parallel::strDist() {
    std::string s;
    for(typename std::map<string,unsigned int>::iterator it=_indexOwner.begin();it!=_indexOwner.end();it++){
        s+=Str("","")+"\n"+it->first+"="+it->second+" "+_ownedIndex[it->first]->strNode();
    }
    return s.substr(0,s.find_last_not_of(" "));
}


void Parallel::addBlock(DerivativeBlock *Block){
    // search grains
    std::stringstream hash;
    hash << Block->cInOut[0]->hash() << Block->cInOut[1]->hash();
    auto it = grainMap.find(hash.str());
    if(it == grainMap.end()){
        grain.push_back(new ParallelGrain());
        grainMap[hash.str()] = grain.back();
        grain.back()->block.push_back(Block);
    }else{
        it->second->block.push_back(Block);
    }
}


void Parallel::addGrain(ParallelGrain *Grain, const string & SendReceive){

    if(SendReceive!="send" and SendReceive!="receive" and SendReceive!="either" and SendReceive!="cpuBalance")
        ABORT("SendReceive allowed values: send,receive,either,cpuBalance. Is: "+SendReceive);

    if(Grain->block.size()==0)return; // empty grain

    // check for existing stripes
    unsigned int iRow=0,iCol=0;
    for(;iRow<cross.size();iRow++)
        if(cross[iRow]->index()==Grain->block[0]->cInOut[1]->idx())break;
    for(;iCol<cross.size();iCol++)
        if(cross[iCol]->index()==Grain->block[0]->cInOut[0]->idx())break;

    // no matching stripe, create new
    if(iCol==cross.size() or iRow==cross.size())cross.push_back(new ParallelCross());

    // send or receive blocks:
    // (1) minimize communication
    // (2) balance load between send and receive branches
    bool send=Grain->leftIndex()->size()<Grain->rightIndex()->size() or
            (Grain->leftIndex()->size()<=Grain->rightIndex()->size() and  cross[iCol]->load(true)<=cross[iRow]->load(true));
    send=SendReceive=="send" or (SendReceive!="receive" and send);

    // override for load balance, disregarding communication
    if(SendReceive=="cpuBalance")send=cross[iCol]->load(true)<cross[iRow]->load(true);

    if(send){
        for(unsigned int k=0;k<Grain->block.size();k++)cross[iCol]->colBlock.push_back(Grain->block[k]);
    }
    else {
        for(unsigned int k=0;k<Grain->block.size();k++)cross[iRow]->rowBlock.push_back(Grain->block[k]);
    }

    if(cross.back()->load()==0.)cross.pop_back(); // newly created cross was not used, remove
}

double Parallel::load() const{
    double l=0;
    if(proc!=0)return l+=proc->load();
    for(unsigned int k=0;k<childSize();k++)l+=child(k)->load();
    return l;
}


static bool smallerCross(const ParallelCross * A, const ParallelCross * B){
    for(size_t k=0;k<A->index()->index().size();k++){
        if(A->index()->index()[k]<B->index()->index()[k])return true;
        if(A->index()->index()[k]>B->index()->index()[k])return false;
    }
    return false;
}

void sortCrosses(std::vector<ParallelCross*> & Cross){
    std::sort(Cross.begin(),Cross.end(),smallerCross);
}


void Parallel::distribute(const string & SendReceive){

    // (re-)arrange ParallelGrain's into ParallelCross's
    for(unsigned int k=0;k<cross.size();k++)delete cross[k];
    cross.clear();
    for(unsigned int k=0;k<grain.size();k++)addGrain(grain[k],SendReceive);

    // get total load
    double loadPerProc=0.;
    for(unsigned int k=0;k<cross.size();k++)loadPerProc+=cross[k]->load(true);
    loadPerProc=loadPerProc/double(MPIwrapper::Size()); // assumes homogenous floating performance

    //for the subinterval
    bool usekGrid = false;

    ParallelLayout lay(cross[0]->index()->root(),"");
    lay.sort(cross);

    double curLoad=0.;
    int curProc=0;

    for(ParallelCross* c: cross)
    {
        if(owner(c->index())==none){
            Parallel::setIndexOwner(c->index(),curProc);
            curLoad+=c->load(true);
            if(curLoad>loadPerProc*(curProc+1)){
                curProc=std::min(curProc+1,MPIwrapper::Size());
            }
#ifdef _DEVELOP_
            int maxProc=curProc;
            int minProc=-curProc;
            MPIwrapper::AllreduceMAX(&maxProc,1);
            MPIwrapper::AllreduceMAX(&minProc,1);
            MPIwrapper::Barrier();
            if(maxProc!=-minProc)PrintOutput::DEVwarning(Sstr+"inconsistent owners(1)"+maxProc+minProc+c->str());
            maxProc= owner(c->nonOwnerIndex());
            minProc=-owner(c->nonOwnerIndex());
            MPIwrapper::AllreduceMAX(&maxProc,1);
            MPIwrapper::AllreduceMAX(&minProc,1);
            if(maxProc!=-minProc)PrintOutput::DEVwarning(Sstr+"inconsistent others(1)"+maxProc+minProc+c->str());
#endif
        } else {
#ifdef _DEVELOP_
            int maxProc= owner(c->index());
            int minProc=-owner(c->index());
            MPIwrapper::AllreduceMAX(&maxProc,1);
            MPIwrapper::AllreduceMAX(&minProc,1);
            if(maxProc!=-minProc){
                PrintOutput::DEVwarning(Sstr+"inconsistent owners(2)"+maxProc+minProc+c->str()+c->index()+c->nonOwnerIndex());
            }
            maxProc= owner(c->nonOwnerIndex());
            minProc=-owner(c->nonOwnerIndex());
            MPIwrapper::AllreduceMAX(&maxProc,1);
            MPIwrapper::AllreduceMAX(&minProc,1);
            if(maxProc!=-minProc)PrintOutput::DEVwarning(Sstr+"inconsistent others(2)"+maxProc+minProc+c->str());
#endif
        }
    }

    // fill up processes
    for(unsigned int k=0;k<cross.size();k++)
    {
        unsigned int nCur(0);
        if(cross[k]->index()->hierarchy().find("k") != string::npos and MPIwrapper::Size() > 1)
            usekGrid = true;
        if(usekGrid){
            setIndexOwner(cross[k]->index(),cross[k]->index()->index()[0]);
            nCur = cross[k]->index()->index()[0];
        }else{
            if(_indexOwner.count(cross[k]->index()->hash())==1){
                // index already has owner
                nCur=owner(cross[k]->index());
            }
            else {
                // find first unfilled process or take lowest filled
                double loadMin=DBL_MAX;
                for(int n=0;n<MPIwrapper::Size();n++){
                    double l=process(n)->load();
                    if(loadMin>l){
                        loadMin=l;
                        nCur=n;
                    }
                    if((n>0 and l<loadPerProc) or l<loadPerProc*0.95){
                        nCur=n;
                        break;
                    }
                }
                setIndexOwner(cross[k]->index(),nCur);
            }
        }

        if(MPIwrapper::isMaster())process(nCur)->addCross(cross[k]);
        int maxProc= owner(cross[k]->index());
        int minProc=-owner(cross[k]->index());
        MPIwrapper::AllreduceMAX(&maxProc,1);
        MPIwrapper::AllreduceMAX(&minProc,1);
        if(maxProc!=-minProc)DEVABORT(Sstr+"inconsistent owners"+maxProc+minProc+cross[k]->index()->strNode()+cross[k]->nonOwnerIndex()->strNode());
        maxProc= owner(cross[k]->nonOwnerIndex());
        minProc=-owner(cross[k]->nonOwnerIndex());
        MPIwrapper::AllreduceMAX(&maxProc,1);
        MPIwrapper::AllreduceMAX(&minProc,1);
        if(maxProc!=-minProc)DEVABORT(Sstr+"inconsistent others"+maxProc+minProc+cross[k]->str()+Parallel::strDist());
    }
    // synchronize index ownership
    syncIndexOwner();

    // on slaves, assign cross to process
    if(not MPIwrapper::isMaster())
        for(size_t k=0;k<cross.size();k++){
            process(_indexOwner[cross[k]->index()->hash()])->addCross(cross[k]);
        }
    for(int k=0;k<MPIwrapper::Size();k++)process(k)->setSendRecv(SendReceive);
}

//OperatorFloor* Parallel::notOnHost=new OperatorDUM();

int Parallel::floorHost(const Index* Idx, const Index* Jdx){
    // get current host (or assign new)
    string iHash=Idx->hash(),jHash=Jdx->hash();
    string hash=iHash+jHash;
    if(_grainHost.count(hash)==0){
        if(_layout)_grainHost[hash]=_layout->floorHost(Jdx);
        if(_layout and _grainHost[hash]==none)_grainHost[hash]=_layout->floorHost(Idx);
        if(_grainHost[hash]==none)_grainHost[hash]=_grainHost.size()%MPIwrapper::Size();
    }
    if(MPIwrapper::Size()>1 and
            Idx->hierarchy().find("ValDer") != string::npos and
            Jdx->hierarchy().find("ValDer") != string::npos
            )_grainHost[hash] = all;
    return _grainHost[hash];
}

int Parallel::moveFloor(OperatorFloor* &Floor, int To){

    int from=int(dynamic_cast<OperatorDUM*>(Floor)==0)?MPIwrapper::Rank():-1;
    MPIwrapper::AllreduceMAX(&from,1);
    if(from<0)return none;
    if(from==To)return from;

    std::vector<int> info(5,0);
    std::vector<complex<double> > buf;
    if(MPIwrapper::Rank()==from){
        if(dynamic_cast<OperatorDUM*>(Floor)!=0)DEVABORT("cannot send dummy floor");
        // already in place - nothing to be done
        if(MPIwrapper::Rank()==To)return from;
        // pack up and send to new host
        Floor->pack(info,buf);
        MPIwrapper::Send(info.data(),info.size(),To,1);
        MPIwrapper::Send(buf.data(),buf.size(),To,2);
        OperatorFloor::replace(Floor); // replace with dummy
    }
    else {
        if(MPIwrapper::Rank()!=To)return from; // not for me
        if(dynamic_cast<OperatorDUM*>(Floor)==0)DEVABORT("cannot overwrite");
        MPIwrapper::Recv(info.data(),info.size(),from,1);
        buf.resize(info[4]);
        MPIwrapper::Recv(buf.data(),buf.size(),from,2);
        OperatorFloor::replace(Floor,info,buf); // replace with received
    }
    return from;
}

OperatorFloor* Parallel::operatorFloor(const Index *iIndex, const Index *jIndex, std::function<OperatorFloor*()> factory){
    int host=floorHost(iIndex,jIndex);
    return host==all or host==MPIwrapper::Rank() ? factory() : new OperatorDUM(0.);
}

void Parallel::addBlocks(std::vector<DerivativeBlock> & Blocks, const std::string & SendReceive){
    // lump together blocks with equal indices into ParallelGrain's
    for(unsigned int k=0;k<Blocks.size();k++)addBlock(&Blocks[k]);
    distribute(SendReceive);
    // create buffers and sender/recipient lists
    for(int n=0;n<MPIwrapper::Size();n++)process(n)->setBuffers(this);
}
unsigned int Parallel::owner(const Index* Idx){
    if(MPIwrapper::Size()==1)return 0;
    if(_indexOwner.count(Idx->hash())==1)return _indexOwner[Idx->hash()];
    return none;
}

/// distribute Coefficients into local
void Parallel::scatter(CoefficientsGlobal *Glob, CoefficientsLocal *Loc, unsigned int From){
    MPIwrapper::ScatterV(Glob->storageData(),Glob->sizes().data(),Loc->storageData(),Loc->size(),From);
}

/// gather Coefficients into Global
void Parallel::gather(CoefficientsGlobal *Glob, CoefficientsLocal *Loc, unsigned int To){
    MPIwrapper::GatherV(Loc->storageData(),Loc->size(),Glob->storageData(),Glob->sizes().data(),To);
}

void Parallel::allGather(CoefficientsGlobal *Glob, CoefficientsLocal *Loc){
    MPIwrapper::AllGatherV(Loc->storageData(),Loc->size(),Glob->storageData(),Glob->sizes().data());
}

/// various tests for parallel operations on Coefficients
void Parallel::test(Coefficients *C){
    MPIwrapper::Barrier();

    Coefficients c(*C);
    CoefficientsGlobal globalC(c.idx());
    CoefficientsGlobal f(c.idx());
    CoefficientsLocal  locC(c.idx());
    c.setToRandom();
    globalC=c;
    if(not (c-=globalC).isZero())ABORT("assignement to global failed");

    c=globalC;
    Parallel::scatter(&globalC,&locC,MPIwrapper::master());
    Parallel::gather(&f,&locC,MPIwrapper::master());
    MPIwrapper::Barrier();

    if(MPIwrapper::isMaster()){
        globalC-=f;
        if(not globalC.isZero())DEVABORT(globalC.str(2)+"\nscatter/gather failed");

        c=f;
        if(not (f-=c).isZero())DEVABORT(globalC.str(2)+"\nback assignement failed");
    }

    MPIwrapper::Barrier();

    globalC.setToRandom();
    // make sure we have the same on all processes
    // why did this ever work w/o this???
    MPIwrapper::AllreduceSUM(globalC.storageData(),globalC.size());
    f=globalC;
    Parallel::keepLocal(&f);
    MPIwrapper::AllreduceSUM(f.storageData(),f.size());
    if(not (globalC-=f).isZero())ABORT("local/sum failed");


    CoefficientsGlobal * viewC=CoefficientsGlobal::view(&c); // data-contiguous global view of c;
    viewC->setToRandom();
    // make sure we have the same on all processes
    MPIwrapper::AllreduceSUM(viewC->storageData(),viewC->size());


    globalC=c;
    CoefficientsLocal *  viewG=CoefficientsLocal::view(&globalC);  // data-contiguous local view on g

    Parallel::gather(viewC,viewG,MPIwrapper::master());

    if(MPIwrapper::isMaster()){
        if(not (c-=globalC).isZero())DEVABORT("global/local views failed");
    }

    Parallel::scatter(&globalC,&locC,MPIwrapper::master());
    if((*viewG-=locC).isZero())
        PrintOutput::DEVmessage(Str("scatter OK"));
    else
        DEVABORT(Str()+MPIwrapper::Rank()+"\n"+viewG->str(1)+"\n"+locC.str(1)+"\nscatter to views failed")
}
void Parallel::gatherAllEigen(const Index *Idx, std::vector<std::complex<double> >& Eval,
                              std::vector<Coefficients *> &EigenVec, std::vector<Coefficients *> &DualVec){

    if(MPIwrapper::Size()==1)return;

    vector<int> siz={(int)Eval.size(),(int)EigenVec.size(),(int)DualVec.size()};
    MPIwrapper::AllreduceSUM(siz.data(),3);
    // end time-consuming
    if(siz[0]==0)return;

    // pack and gather
    ParallelBuffer bufE,bufV,bufD;
    bufE.add(Eval);bufE.allGather();
    for(size_t k=0;k<EigenVec.size();k++)bufV.add(EigenVec[k]);
    for(size_t k=0;k<DualVec.size();k++) bufD.add( DualVec[k]);
    bufV.allGather();
    bufD.allGather();

    // extend storage
    Eval.resize(siz[0]);
    for(int k=EigenVec.size();k<siz[1];k++)EigenVec.push_back(new Coefficients(Idx));
    for(int k= DualVec.size();k<siz[2];k++) DualVec.push_back(new Coefficients(Idx));

    // overwrite by gathered data
    bufE.extract(0,Eval);
    for(int k=0,pos=0;k<siz[1];k++)pos=bufV.extract(pos,EigenVec[k]);
    for(int k=0,pos=0;k<siz[2];k++)pos=bufD.extract(pos, DualVec[k]);
}

void Parallel::unsetIndexOwner(const Index *Idx){
    auto pOwn=_indexOwner.find(Idx->hash());
    if(pOwn!=_indexOwner.end())_indexOwner.erase(pOwn);
    auto qOwn=_ownedIndex.find(Idx->hash());
    if(qOwn!=_ownedIndex.end())_ownedIndex.erase(qOwn);
}

void Parallel::setIndexOwner(const Index *Idx, int Proc){
    if(MPIwrapper::Size()>1 and owner(Idx)!=Parallel::none)DEVABORT("Index has owner, cannot be set");
    _indexOwner[Idx->hash()]=Proc;
    _ownedIndex[Idx->hash()]=Idx;
}

void Parallel::syncIndexOwner(){

    vector<int> own;
    for(map<std::string,unsigned int>::iterator it=_indexOwner.begin();it!=_indexOwner.end();it++)
        own.push_back(it->second);

    own.push_back(_indexOwner.size());

    MPIwrapper::Bcast(own.data(),own.size(),MPIwrapper::master());

    if(int(_indexOwner.size())!=own.back())ABORT(Str("failed basic consistency requirement: index size differs, on master=")
                                                 +own.back()+"on thread="+_indexOwner.size());
    int k=0;
    for(map<std::string,unsigned int>::iterator it=_indexOwner.begin();it!=_indexOwner.end();it++,k++)
        it->second=own[k];
}


void Parallel::timer(){
    string file=ReadInput::main.output()+"timer";
    ofstream out(file.c_str());
    Timer::write(out);
}
