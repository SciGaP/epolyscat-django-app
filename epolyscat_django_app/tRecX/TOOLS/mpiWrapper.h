// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __MPIWRAPPER_H__
#define __MPIWRAPPER_H__

#include <vector>
#include <string>
#include <complex>
#include <iostream>
#include <fstream>
#include "limits.h"
#include "str.h"

#ifndef _NOMPI_
#define __USE_MPI__
#define __NO_SHM__
#else
typedef int MPI_Request;
typedef int MPI_Comm;
typedef int MPI_Group;
typedef int MPI_Group;
#endif //_NOMPI_

#ifndef _NOMPI_
#include <mpi.h>
#ifndef __NO_SHM__
//    #include <boost/interprocess/managed_shared_memory.hpp>
//    #include <boost/interprocess/containers/vector.hpp>
//    #include <boost/interprocess/allocators/allocator.hpp>
//    typedef boost::interprocess::allocator<std::complex<double>, boost::interprocess::managed_shared_memory::segment_manager> ShmemAllocator;
//    typedef boost::interprocess::vector<std::complex<double>, ShmemAllocator> MyVector;
#endif
#endif

#define MPout MPIwrapper::Out()
#define Pout MPIwrapper::strRank; ///< print with Rank
#define MASTEROUT if(MPIwrapper::isMaster())MPIwrapper::Out()
#define Mout Str("","")+"<"+MPIwrapper::Rank(MPIwrapper::worldCommunicator())+">"+Str::sep(" ")


namespace MPIwrapper
{
MPI_Comm communicator();
MPI_Comm worldCommunicator();
MPI_Comm setCommunicator(const std::vector<int> ProcRanks, MPI_Comm ParentCommunicator=MPIwrapper::communicator());
MPI_Comm setCommunicator(MPI_Comm ParentCommunicator);

const int undefined=INT_MAX;

void print(std::string String);

class Buffer{
public:
    virtual ~Buffer(){}
    Buffer():req(0){}
    MPI_Request req;
    std::vector<std::complex<double> > val;
    unsigned int size() const {return val.size();}
};


bool isMPI();
void Init(int argc, char* argv[]);
void Finalize();
void Exit(int i);
void Abort(int i);
int Rank(MPI_Comm Comm=MPIwrapper::communicator());
int Size(MPI_Comm Comm=MPIwrapper::communicator());

std::string Get_processor_name();
void Barrier(MPI_Comm Comm=MPIwrapper::communicator());
void Barrier_slowButHard();
void Send(int & x, int dest, int tag);
void Recv(int & x, int dest, int tag);
void Send(int * x, int size, int dest, int Tag=undefined);
int Recv(int * x, int size, int dest, int Tag=undefined);
void Send(char * x, int size, int dest, int Tag=undefined);
void Recv(char * x, int size, int dest, int Tag=undefined);
void Send(double & x, int dest, int tag);
void Recv(double & x, int dest, int tag);
void Send(std::string & Str, int dest);
int Recv(std::string & Str, int dest);

void GatherV(   std::complex<double>*Loc, int LSize, std::complex<double>* Glob, int * GSizes,  int From=0);
void AllGatherV(std::complex<double> *Loc, int LSize, std::complex<double> *Glob, int *GSizes);
void ScatterV(std::complex<double>* Glob, int * GSizes, std::complex<double>*Loc, int LSize, int From=0);

// raw Bcast
void Bcast(char* Data, int Size, int From);
void Bcast(int* Data, int Size, int From);
void Bcast(double* Data, int Size, int From);
void Bcast(std::complex<double>* Data, int Size, int From);

// raw initiate send/recv
void ISend(std::complex<double> * Data, int Size, int To,   MPI_Request & Req, int Tag=MPIwrapper::undefined);
void IRecv(std::complex<double> * Data, int Size, int From, MPI_Request & Req, int Tag=MPIwrapper::undefined);

// raw recv
void Recv(std::complex<double> * Data, int Size, int From, int Tag=MPIwrapper::undefined);
void Send(std::complex<double> *Data, int Size, int To, int Tag=MPIwrapper::undefined);

// raw reduce (aliasing)
void AllreduceSUM(std::complex<double> * Data, int Siz);
void AllreduceSUM(int * Data, int Siz);
void AllreduceSUM(int & Data);
void AllreduceSUM(double * Data, int Siz);
void AllreduceMAX(int * Data, int Siz);
void AllreduceMAX(long &Data);
void AllreduceMAX(long long & Data);
void AllreduceMAX(double * Data, int Siz);

// secondary convenience functions
void Bcast(std::string & Data, int From);
void Bcast(std::vector<std::string> &Data, int From);

void Waitall(const std::vector<Buffer> & Buf);


void AllreduceMAX(double & send, double & recv);
void AllreduceSUM(std::complex<double> & send, std::complex<double> & recv);

int  master();
int  anySource();
bool isMaster(MPI_Comm Comm=MPIwrapper::communicator());

class Out {
    bool start;
public:
    Out():start(true){}
    template<class T>
    std::ostream& operator<<(const T& x){
        if(start){
            std::cout<<"<"<<MPIwrapper::Rank()<<"> ";
        }
        start=false;
        std::cout<<x;
        return std::cout;
    }
    Out& operator<<(std::ostream& (*F)(std::ostream&)){
        F(std::cout);
        return *this;
    }
};
void printAll(std::string S);

class File : public std::ofstream {
public:
    ~File(){}
    File(const std::string Name);
};

void ISend(std::vector<std::complex<double> > & coef, int dest, int tag, MPI_Request *r);
void IRecv(std::vector<std::complex<double> > & coef, int dest, int tag, MPI_Request *r);
void Wait(MPI_Request * r);
void Waitall(std::vector<MPI_Request> & array_of_requests);
void ISend(int & x, int dest, int tag, MPI_Request *r);
void IRecv(int & x, int dest, int tag, MPI_Request *r);
void ISend(std::vector<int> & x, int dest, int tag, MPI_Request *r);
void IRecv(std::vector<int> & x, int dest, int tag, MPI_Request *r);

#ifndef _NOMPI_
#ifndef __NO_SHM__
class SharedMemoryXcdVector{
public:
    SharedMemoryXcdVector(int shared_memory_name, int size, std::string name="", int id_to_create=0, bool cout_info=true);
    ~SharedMemoryXcdVector();

    std::complex<double> * data();
    unsigned int size();

private:
    MyVector * x;
    std::string shared_memory_name;
    boost::interprocess::managed_shared_memory* segment;
    //           MyShared* segment;
};
#endif
#endif
}

#endif // MPIWRAPPER_H
