// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "parallelContinuity.h"

#include "timer.h"

#include "basisIntegrable.h"
#include "coefficientsLocal.h"
#include "index.h"
#include "parallel.h"
#include "printOutput.h"
#include "mpiWrapper.h"
using namespace std;

map<string, vector<map<string,vector<ParallelContinuity::Marg> > > >ParallelContinuity::neighbor;

ParallelContinuity::ParallelContinuity(Coefficients * C, int AtBoundary)
{
    if(C==0)ABORT("cannot not construct parallel contininuity for zero coefficient pointer");
    if(C->size()==0)return; // empty coefficient - no need to
    string hashC=C->idx()->hash()+"_"+tools::str(AtBoundary);
    if(neighbor.count(hashC)==0)
        setIndex(hashC,C->idx(),0,AtBoundary);
    pMarg.resize(neighbor[hashC].size());
    locCur.resize(neighbor[hashC].size());
    locNei.resize(neighbor[hashC].size());
    sendBuf.resize(neighbor[hashC].size());
    recvBuf.resize(neighbor[hashC].size());
    for(unsigned int d=0;d<neighbor[hashC].size();d++){
        // hybrid--------
        if(C->idx()->isHybrid()){
            std::vector<std::string>subs=tools::splitString(C->idx()->axisSubset(),'&');
            if(subs.size()!=2)DEVABORT("Idx.isHybrid, but axisSubset="+C->idx()->axisSubset());
            Coefficients* cRoot=C;
            for(size_t k=0;k<cRoot->childSize();k++)
                if(C->child(k)->idx()->axisSubset()==subs[1])
                    C=cRoot->child(k);
        }
        //--------------------
        // continuity for each direction
        Coefficients * cCur=C->firstLeaf();
        pMarg[d].resize(MPIwrapper::Size());
        // change
        for(unsigned int i = 0; i < pMarg[d].size(); i++)
            pMarg[d][i].resize(2); // one for lower margin, the other for upper margin
        for(;cCur!=0;cCur=cCur->nodeRight()){
            //NOTE: sequence matters, from lower to upper
            unsigned int oCur=Parallel::owner(cCur->idx());
            for(unsigned int k=0;k<neighbor[hashC][d][cCur->idx()->hash()].size();k++)
            {
                // information about margin layout of current Coefficient leaf
                Marg marg=neighbor[hashC][d][cCur->idx()->hash()][k];
                unsigned int oNei=Parallel::owner(marg.fNei);
                unsigned int pcur=marg.curPosFloor();
                unsigned int pnei=marg.neiPosFloor();
                if(oNei!=oCur){
                    unsigned int lowUp=marg.curPosIndex() < marg.neiPosIndex();
                    for(int k=0;k<marg.curSize();k++)pMarg[d][oNei][lowUp].push_back(cCur->floorData()+pcur+k);
                }
                else {
                    if(marg.fNei->posIndex()>cCur->idx()->posIndex()){
                        // margins with neighbors on local node
                        // upper and lower margin are done in one shot
                        //                        Coefficients * cNei=C->root()->retrieve(marg.fNei); // get neighbor Coefficients
                        Coefficients * cNei=C->retrieve(marg.fNei); // get neighbor Coefficients
                        for(int k=0;k<marg.curSize();k++){
                            locCur[d].push_back(cCur->floorData()+pcur+k);
                            locNei[d].push_back(cNei->floorData()+pnei+k);
                        }
                    }
                }
            }
        }
        // send and receive buffers
        sendBuf[d].resize(MPIwrapper::Size());
        recvBuf[d].resize(MPIwrapper::Size());

        for(unsigned int n=0;n<sendBuf[d].size();n++){
            sendBuf[d][n].resize(2);
            recvBuf[d][n].resize(2);
            for(unsigned int i = 0; i < sendBuf[d][n].size(); i++){
                sendBuf[d][n][i].resize(pMarg[d][n][i].size());
                recvBuf[d][n][i].resize(pMarg[d][n][i].size());
            }
        }
    }
}

// get neighbor indices
void ParallelContinuity::setIndex(string HashC, const Index *I, unsigned int Dimension, int AtBoundary){

    if(I->continuity()!=Index::npos){
        if(neighbor[HashC].size()<=Dimension)neighbor[HashC].resize(Dimension+1);
        // get lower and upper Marg's on all boundaries between branches
        int begK=1,endK=I->childSize();
        if(AtBoundary!=-1){begK=AtBoundary;endK=begK+1;}
        for(int k=begK;k<endK;k++){
            //NOTE: sequence matters: first do lower margin (i.e. from k down)
            addNeighbor(HashC,I->child(k),  I->child(k-1),Dimension,I->continuity());
            addNeighbor(HashC,I->child(k-1),I->child(k),  Dimension,I->continuity());
        }
        if(AtBoundary!=-1)return; // only a single split at first FEM level encountered
        Dimension++;
    }

    // recursively get neighbors on further levels
        for(size_t k=0;k<I->childSize();k++){
            if(not I->child(k)->isLeaf())
                setIndex(HashC,I->child(k),Dimension,AtBoundary);
        }
}

void ParallelContinuity::addNeighbor(std::string HashC, Index *Cur, Index *Nei, unsigned int Dimension, unsigned int Level){


    if(Cur->depth()==Level){
        // on the continuity level

        // get the floor indices
        const Index * fCur=Cur,*fNei=Nei;
        while(not fCur->hasFloor())fCur=fCur->parent();
        while(not fNei->hasFloor())fNei=fNei->parent();
        unsigned int mCur=Cur->basis()->integrable()->upperMargin();
        unsigned int mNei=Nei->basis()->integrable()->lowerMargin();
        if(Cur->posIndex()>Nei->posIndex()){
            mCur=Cur->basis()->integrable()->lowerMargin();
            mNei=Nei->basis()->integrable()->upperMargin();
        }
        // insert  indices in floor neighbor and floor index into table
        //        neighbor[HashC][Dimension][fCur->hash()].push_back(Marg(Cur->child(mCur),Nei->child(mNei),fNei));
        neighbor[HashC][Dimension][fCur->hash()].push_back(Marg(Cur,mCur,Nei,mNei,fNei));
    }
    else {
        for(unsigned int k=0;k<Cur->childSize();k++)
            addNeighbor(HashC,Cur->child(k),Nei->child(k),Dimension,Level);
    }
}

TIMER(parCont,)
TIMER(parContA,)
TIMER(parContB,)
TIMER(parContC,)
void ParallelContinuity::apply(Coefficients *C,double Scal){

    STARTDEBUG(parCont);
    // need to do directions sequentially
    for(unsigned int d=0;d<sendBuf.size();d++){

        // initiate send and recveive
        vector<vector<MPI_Request> > sendReq(sendBuf[d].size());
        vector<vector<MPI_Request> > recvReq(recvBuf[d].size());
        for(unsigned int n=0;n<sendBuf[d].size();n++){
            sendReq[n].resize(2);
            recvReq[n].resize(2);
            // place margin values into send buffer
            // split the boundaries
            for(unsigned int i = 0; i < 2; i++){
                for(unsigned int k=0;k<sendBuf[d][n][i].size();k++)sendBuf[d][n][i][k]=*pMarg[d][n][i][k];
                if(sendBuf[d][n][i].size()>0)MPIwrapper::ISend(sendBuf[d][n][i].data(),sendBuf[d][n][i].size(),n,sendReq[n][i],sendBuf[d][n][i].size());
                if(recvBuf[d][n][i].size()>0)MPIwrapper::IRecv(recvBuf[d][n][i].data(),recvBuf[d][n][i].size(),n,recvReq[n][i],sendBuf[d][n][i].size());
            }

        }
        for(size_t k=0;k<locCur[d].size();k++){
            *locNei[d][k]=*locCur[d][k]=0.5*(*locNei[d][k]+*locCur[d][k])*Scal;
        }

        // wait for receive and average with remote margins
        for(int n=0;n<MPIwrapper::Size();n++){
            //NOTE: we must ALWAYS match non-blocking operations with wait's, else a memory leak
            for(unsigned int i = 0; i < 2; i++){
                if(sendBuf[d][n][i].size()>0)MPIwrapper::Wait(&sendReq[n][i]);
                if(recvBuf[d][n][i].size()>0)MPIwrapper::Wait(&recvReq[n][i]);
                for(unsigned int k=0;k<recvBuf[d][n][i].size();k++)*pMarg[d][n][i][k]=0.5*(*pMarg[d][n][i][k]+recvBuf[d][n][i][k])*Scal;
            }

        }
    }
    STOPDEBUG(parCont);
}

void ParallelContinuity::margin(ParallelContinuity * Marg) const
{
    for(unsigned int d=0;d<locCur.size();d++){

        // local margins
        for(unsigned int k=0;k<locCur[d].size();k++)
            *Marg->locCur[d][k]=*locCur[d][k];

        // shared margins
        for(unsigned int n=0;n<MPIwrapper::Size();n++){
            for(unsigned int i = 0; i < 2; i++)
                for(unsigned int k=0;k<pMarg[d][n][i].size();k++)
                    *Marg->pMarg[d][n][i][k]=*pMarg[d][n][i][k];
        }
    }
}

void ParallelContinuity::setMargin(complex<double> Val)
{
    for(unsigned int d=0;d<locCur.size();d++){

        // local margins
        for(unsigned int k=0;k<locCur[d].size();k++)
            *locNei[d][k]=*locCur[d][k]=Val;

        // shared margins
        for(unsigned int n=0;n<MPIwrapper::Size();n++){
            for(unsigned int i = 0; i < 2; i++)
                for(unsigned int k=0;k<pMarg[d][n][i].size();k++)
                    *pMarg[d][n][i][k]=Val;
        }
    }
}
void ParallelContinuity::halfDiffMargin(ParallelContinuity * Marg) const
{

    // need to do directions sequentially
    for(unsigned int d=0;d<Marg->sendBuf.size();d++){

        // initiate send and recveive
        vector<vector<MPI_Request> > sendReq(Marg->sendBuf[d].size());
        vector<vector<MPI_Request> > recvReq(Marg->recvBuf[d].size());
        for(unsigned int n=0;n<Marg->sendBuf[d].size();n++){
            // place margin values into send buffer
            sendReq[n].resize(pMarg[d][n].size());
            recvReq[n].resize(pMarg[d][n].size());
            for(unsigned int i = 0; i < pMarg[d][n].size(); i++){
                for(unsigned int k=0;k<sendBuf[d][n][i].size();k++)Marg->sendBuf[d][n][i][k]=*pMarg[d][n][i][k];
                if(Marg->sendBuf[d][n][i].size()>0)MPIwrapper::ISend(Marg->sendBuf[d][n][i].data(),Marg->sendBuf[d][n][i].size(),n,sendReq[n][i]);
                if(Marg->recvBuf[d][n][i].size()>0)MPIwrapper::IRecv(Marg->recvBuf[d][n][i].data(),Marg->recvBuf[d][n][i].size(),n,recvReq[n][i]);
            }
        }

        // diff local margins
        for(unsigned int k=0;k<locCur[d].size();k++){
            *Marg->locCur[d][k]=0.5*(*locCur[d][k]-*locNei[d][k]);
            *Marg->locNei[d][k]=-*Marg->locCur[d][k];
        }

        // wait for receive and average with remote margins
        for(unsigned int n=0;n<MPIwrapper::Size();n++){
            //NOTE: we must ALWAYS match non-blocking operations with wait's, else a memory leak
            for(unsigned int i = 0; i < Marg->sendBuf[d][n].size(); i++){
                if(Marg->sendBuf[d][n][i].size()>0)MPIwrapper::Wait(&sendReq[n][i]);
                if(Marg->recvBuf[d][n][i].size()>0)MPIwrapper::Wait(&recvReq[n][i]);
                for(unsigned int k=0;k<Marg->recvBuf[d][n][i].size();k++)
                    *Marg->pMarg[d][n][i][k]=0.5*(*pMarg[d][n][i][k]-Marg->recvBuf[d][n][i][k]);
            }

        }
    }
}
ParallelContinuity::Marg::Marg(const Index * PCur,int NCur, const Index * PNei, int NNei, const Index * FNei)
    :_cur(0),_nei(0),fNei(FNei)
{
    if(not PCur->isBottom()){
        pCur=PCur->child(NCur)->posInFloor();
        pNei=PNei->child(NNei)->posInFloor();
        iCur=PCur->child(NCur)->posIndex();
        iNei=PNei->child(NNei)->posIndex();
        _size=PCur->child(NCur)->sizeStored();
    }
    else {
        pCur=const_cast<Index*>(PCur)->posInFloor()+NCur;
        pNei=const_cast<Index*>(PNei)->posInFloor()+NNei;
        iCur=const_cast<Index*>(PCur)->posIndex()+NCur;
        iNei=const_cast<Index*>(PNei)->posIndex()+NNei;
        _size=1;
    }
}


ParallelContinuity::Marg::Marg(const Index * Cur,const Index * Nei, const Index * FNei)
    :_cur(Cur),_nei(Nei),fNei(FNei){
    pCur=Cur->posInFloor();
    pNei=Nei->posInFloor();
    iCur=Cur->posIndex();
    iNei=Nei->posIndex();
    _size=Cur->sizeStored();
}
int ParallelContinuity::Marg::curPosFloor(){return pCur;}
int ParallelContinuity::Marg::neiPosFloor(){return pNei;}
int ParallelContinuity::Marg::curPosIndex(){return iCur;}
int ParallelContinuity::Marg::neiPosIndex(){return iNei;}
int ParallelContinuity::Marg::curSize(){return _size;}
