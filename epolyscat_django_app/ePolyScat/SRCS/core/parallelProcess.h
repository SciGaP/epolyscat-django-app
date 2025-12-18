// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARALLELPROCESS_H
#define PARALLELPROCESS_H

#include <vector>
#include <complex>
#include <deque>
#include <map>
#include "parallelGrain.h"

class Index;
class Parallel;
//class ParallelGrain;
class ParallelCommunicate;
class ParallelContinuity;
class ParallelPacket;
class ParallelCross;

#include "mpiWrapper.h"

/** \ingroup Parallelization */
/// parallel structures on one process
class ParallelProcess
{
    friend class Parallel;
    friend class DerivativeFlat;
    friend class ParallelOperator;

    unsigned int numb;            // unique node number (usually == node rank)
    bool _used; // indicates non-zero process

    std::vector<Coefficients*> _outFloors; ///< pointers lhs floors touched by present process

    std::vector<const DerivativeBlock*> temp;                  ///< temporary block pointer storage
    std::vector<const DerivativeBlock*> localB;                ///< blocks with local in- and output
    std::deque<std::vector<const DerivativeBlock*> > recvB;   ///< blocks that receive input
    std::deque<std::vector<const DerivativeBlock*> > sendB;   ///< blocks that send output

    void check(); // various checks afte setup of ParallelProcess is complete

    std::vector<std::vector<Coefficients*> > recvC; ///< before apply, get these coefficients
    std::vector<std::vector<Coefficients*> > sendC; ///< after  apply, send these coeffiecients

    std::vector<MPIwrapper::Buffer> recvBuf; ///< before apply, recvB recvBuf[k] from process(k)
    std::vector<MPIwrapper::Buffer> sendBuf; ///< after apply, send sendBuf[k] to process(k)

    void sendTo(unsigned int Recipient, std::vector<Coefficients*> & C, MPIwrapper::Buffer &Buf);
    void recvFrom(unsigned int Sender,  std::vector<Coefficients*> & C, MPIwrapper::Buffer & Buf);
    void  addFrom(unsigned int Sender,  std::vector<Coefficients*> & C, MPIwrapper::Buffer & Buf);
public:
    static bool debugTimer; ///< for activating selected timers;
    ~ParallelProcess();
    void addCross(ParallelCross* Cross); ///< sort DerivativeBlock's of Cross into ParallelProcess::send,recvB,local
    void setSendRecv(const std::string &Assign="local"); ///< local/send/recv blocks
    ParallelProcess(unsigned int Nproc, unsigned int Numb);
    unsigned int number() const {return numb;}
    bool unused() const;
    double load() const;
    void setBuffers(const Parallel * Par);
    std::string str() const;
};

#endif // PARALLELPROCESS_H
