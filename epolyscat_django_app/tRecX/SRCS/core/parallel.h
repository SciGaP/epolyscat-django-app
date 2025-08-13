// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARALLEL_H
#define PARALLEL_H

#include <complex>
#include <vector>
#include <map>
#include <memory>
#include "tree.h"
//#include "derivativeFlat.h"
#include "parallelCross.h"
#include "parallelLayout.h"
#include "coefficients.h"
#include "coefficientsFloor.h"
#include "coefficientsLocal.h"
#include "index.h"
#include "threads.h"

class OperatorAbstract;
//class Operator;
class OperatorFloor;
class Discretization;

class ParallelProcess;
class ParallelCross;
class ParallelGrain;
class ParallelLayout;

class MPI_request;


/** @defgroup Parallelization
 * \brief distribution, load-balancing, communication, hardware-structure, etc.
 *  @{
*/

///@brief Central structure for parallelization
///
/// Basic principles
///
/// - the parallel layout is organized into a Tree<Parallel>
///   whose lowest level contains ParallelProcess's
/// <br>the Tree levels above can be used to encode hardware layout
/// <br>at the moment we distiguish
/// <br>"Nodes" ...asssumed to be slowly linked (e.g. Infiniband)
/// <br>"Boards"...separate boards on the same node
/// <br>"CPUs"  ...separate CPUs on the same board
/// <br>"Cores" ...separate cores within a CPU on the same board
/// <br>(this can be taken into account in Parallel::distribute(),
///    but is not exploited at present)
///
/// - for any given Index, each "floor" will owned by exactly one ParallelProcess
/// <br>a floor is a subtree whose data will be stored contiguously
///   and it is uniquely determined by its position in the Index tree
///
/// - for parallelization, the sequence of floors may be redefined by permuting
///   their multi-index, e.g. such that finite element indices come first
///   for a given Parallel object, this is set by Parallel::setSort
/// <br>sorting of the multi-indices (i0,i1,..,iN) is row-wise, i.e. i0 runs slowest
///
/// - any operator is split into DerivativeBlock's
///   (a DerivativeBlock is a map between floors)
///
/// - there is one ParallelCross each floor: this is a set of blocks where either left or right Index is the given floor
/// <br> a ParallelCross is the smallest unit that can be assigned to a ParallelProcess
/// <br> ParallelCross's are set up to all have approximately equal ParallelCross::load()
///
/// - Parallel::distribute() assigns ParallelCross's to ParallelProcess'es using ParallelCross::load()
/// <br>(see Parallel::distribute() for the algorithm)
///
/// - during Parallel::distribute(), Index floor ownership is assigned to the ParallelProcess
/// <br>(see Parallel::assignFloorOwner() for the algorithm)
///
/// - going through through all ParallelProcess'es, each process sets up the sends and receives it needs
///

class Parallel:public Tree<Parallel>{
    friend class ParallelProcess;

    static unsigned int sizeNode;
    static unsigned int sizeBoard;
    static unsigned int sizeCPU;

    void construct(unsigned int NProc, unsigned int Level);

    void addBlock(DerivativeBlock *Block); ///< add a single block to grain

    void distribute(const std::string & SendReceive); ///< assign ParallelCross's  to process
    void reDistribute(); ///< re-distribute the floor data as needed

    std::map<std::string, ParallelGrain*> grainMap;
    //std::vector<ParallelGrain*> grain; /// blocks with equal left and right indices lumped together
    std::vector<ParallelCross*> cross; /// sets of blocks sharing one single index, either right or left
    static std::unique_ptr<ParallelLayout> _layout;
public: //HACK
    std::vector<ParallelGrain*> grain; /// blocks with equal left and right indices lumped together
    std::vector<ParallelProcess*> _process; /// fast access to processes

    /// unique registry for Index's ownership
    static std::map<std::string,unsigned int> _indexOwner;
    static void unsetIndexOwner(const Index* Idx);
    static void setIndexOwner(const Index* Idx, int Proc);
private: //HACK end
    static std::map<std::string,unsigned int> _operatorHost;
    static std::map<std::string,unsigned int> _grainHost;


    /// add a ParallelGrain to suitable ParallelCross
    void addGrain(ParallelGrain *Grain, const std::string &SendReceive);
    void applySubset(const std::vector<const DerivativeBlock *> &Block,std::complex<double>Alfa) const;

    void syncIndexOwner();

public:
    /// current master for OperatorFloor calculation
//    static  std::map<const OperatorFloor*,unsigned int> _floorHost;
    static bool debugTimer; ///< for activating selected timers;


    static void clear(){_indexOwner.clear();_operatorHost.clear();_grainHost.clear();};
    static const unsigned int none=INT_MAX;
    static const unsigned int all=INT_MAX-1;
    static const unsigned int thread=INT_MAX-2;

    ParallelProcess * proc;
    static unsigned int owner(const Index* Idx);//{if(_indexOwner.count(Idx->hash())==1)return _indexOwner[Idx->hash()];return none;}
    static void scatter(CoefficientsGlobal* Glob, CoefficientsLocal *Loc, unsigned int From);
    static void gather(CoefficientsGlobal *Glob, CoefficientsLocal *Loc, unsigned int To);
    static void allGather(CoefficientsGlobal *Glob, CoefficientsLocal *Loc);

    static OperatorFloor* operatorFloor(const Index* iIndex, const Index* jIndex, std::function<OperatorFloor*()> factory);

    static int moveFloor(OperatorFloor *&Floor, int To); ///< move a floor from its current process to new
    static int floorHost(const Index* Idx, const Index* Jdx);
    static void gatherAllEigen(const Index* Idx, std::vector<std::complex<double> > &Eval, std::vector<Coefficients*> &EigenVec, std::vector<Coefficients*> &DualVec);
    static void bCast(OperatorFloor *&Of); ///< broadcast OperatorFloor to all processes

    static void keepLocal(Coefficients *C); ///< for testing: zero all but local data
    static void test(Coefficients* C);

    static void setSort(const Discretization* Disc); ///< using axes, define re-sorting of indices for parallelization
    static std::map<const Index*,std::vector<unsigned int> > indexSort; ///< how a given tree-index is to be sorted
    ~Parallel();//{for(unsigned int k=0;k<cross.size();k++)delete cross[k];}

    static void timer();

    Parallel(ParallelProcess* Proc=0):proc(Proc){}
    Parallel(unsigned int NProc,unsigned int Level):proc(0){construct(NProc,Level);}
    inline ParallelProcess *process(unsigned int N) const {return _process[N];} ///< pointer to node process N (not fast)
    double load() const; ///< load on layout

    /// Load-balanced distribution of Blocks
    void addBlocks(std::vector<DerivativeBlock> &Blocks,
                   const std::string &SendReceive /** for non-local apply, select to first "send", "receive", "either" */);

    void setToZeroLHS(); ///< all lhs's = 0
    // non-const ONLY for development (local re-assignement of proc
    void apply(const std::complex<double> & Alfa); ///< apply operator in parallel layout

    using Tree::str; // this shuts up compiler warning about str() w/o arguments hiding str(...) from Tree
    std::string str() const; ///< print the layout
    static std::string strDist(); ///< print detailed distribution


};

/** @} */ // end group Parallelization

#endif // PARALLEL_H
