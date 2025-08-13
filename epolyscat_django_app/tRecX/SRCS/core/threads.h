// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef THREADS_H
#define THREADS_H

#include "mpiWrapper.h"
#include "coefficients.h"

// do we need a class??
#include <vector>
#include <complex>
#include <memory>

class Index;

#include "basisSub.h"
///@brief Used exclusively as a special BasisSub for threads
class BasisThread: public BasisSub {
public:
    BasisThread(std::string Def):BasisSub("Subset: "+Def.substr(8)){}
    static std::string strDefinition(const BasisAbstract* Bas, std::vector<int> Subset){
        return "Thread: "+BasisSub::strDefinition(Bas,Subset).substr(10);
    }
    std::string str(int Level) const{return "thread"+BasisSub::str(Level).substr(6);}
};


///@ brief Handles threaded Index's and Coefficients
///
/// a given index is split at the floor level
class Threads
{
    int _rows; ///< number of rows in floor (we have row-wise storage: increment if split column-wise)
    std::shared_ptr<const Index> _threadedI; ///< threaded index
    std::shared_ptr<const Index> _joinedI; ///< joint index
    std::shared_ptr<Coefficients> _C; /// joint coefficients (i.e. w/o Threads on top)
    std::shared_ptr<Coefficients> _T; /// threaded coefficients (i.e. with Threads on top)

    const Index* join(std::vector<const Index*> Idx);
    static void join(Coefficients* J, std::vector<const Coefficients*> & VF, int Cols);
    static void scatter(const Coefficients* J, std::vector<Coefficients*> & VF,int Rows);
    static bool set(const Coefficients *C); ///< create Threads object on world master (false if failure)
public:
    ~Threads(){}
    Threads(){}
    Threads(const Index* Idx);
    Threads(const Coefficients* C);

    static void setup(MPI_Comm Communicator);///< create all communicators, must be called by all in MPI_COMM_WORLD
    static MPI_Comm all();     // all Threads
    static MPI_Comm single();  // single Thread (= world Rank)
    static MPI_Comm current(); // current Thread

    static const Index* join(const Index* Idx);
    ///@brief merge Threads into a single Coefficient with top level removed
    ///
    /// storage for C is created upon first call
    static Coefficients *join(Coefficients &C);

    ///@brief new index forked at axis name - top level will be "Threads"+AxisName
    ///
    /// AxisName must be in floor (this is why it is difficult)
    /// <br> each thread contains a share of the AxisName level
    static Index * fork(const Index* Idx, std::string AxisName);

    ///@brief scatters a joint Coefficient C into Scattered
    static void scatter(const Coefficients* C, Coefficients &Scattered);

    /// return thread(rank) as new Index
    static Index * detach(const Index* Idx);

    ///@brief index is branch of Threads
    static bool isThread(const Index* Idx);

    ///@brief master of the Threads
    static bool isMaster();
    static int rank();

    static double max(double Val);
    static double sum(double Val);
    static std::complex<double> sum(std::complex<double> Val);

    static void kill(std::string Message);

};
#endif // THREADS_H
