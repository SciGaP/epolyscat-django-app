// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "mpiWrapper.h"
#include <string>
#include <cstring>

#include "tools.h"
#include "str.h"

#ifdef _WIN32
// windows stuff
#include <io.h>
#else
#include <unistd.h> //for Vinays computer
#endif

//resolve forward declarations
#ifdef _USE_BOOST_
#include <boost/lexical_cast.hpp>
#endif

#ifndef _NOMPI_
#include <mpi.h>
#endif

#ifndef _TIMER_OFF_
#include "timer.h"
static TimerRecursive MPITimer("MPI","MPI",__FILE__);
#define MPISTART MPITimer.start();
#define MPISTOP  MPITimer.stopTimer();
#else
#define MPISTART
#define MPISTOP
#endif

using namespace std;
using namespace tools;

static int MPIwrapper_master=0;
MPI_Comm MPIwrapper::worldCommunicator() {
#ifndef _NOMPI_
    return MPI_COMM_WORLD;
#else
    return 0;
#endif
}
MPI_Comm currentCommunicator=MPIwrapper::worldCommunicator();

MPI_Comm MPIwrapper::communicator() {if(currentCommunicator)return currentCommunicator; return worldCommunicator(); }

MPI_Comm MPIwrapper::setCommunicator(MPI_Comm Communicator){
    currentCommunicator=Communicator;
    return currentCommunicator;
}
MPI_Comm MPIwrapper::setCommunicator(const std::vector<int> ProcRanks, MPI_Comm ParentCommunicator){
#ifndef _NOMPI_
    int err;
    MPI_Group groupParent;
    MPI_Group groupComm;
    MPI_Comm_group(ParentCommunicator, &groupParent);
    std::vector<int>procRanks(ProcRanks);
    err=MPI_Group_incl(groupParent, procRanks.size(), procRanks.data(), &groupComm); // create the new group
    if(err==MPI_ERR_ARG)DEVABORT("Invalid arguments");
    if(err==MPI_ERR_GROUP)DEVABORT("Null or invalid group");
    if(err==MPI_ERR_INTERN)DEVABORT("Internal");
    if(err==MPI_ERR_RANK)DEVABORT("Invalid source or destination rank");
    if(err!=MPI_SUCCESS)DEVABORT("Unclassified error");

    err=MPI_Comm_create(ParentCommunicator,groupComm,&currentCommunicator);
    if(err==MPI_ERR_COMM)DEVABORT("Invalid communicator");
    if(err==MPI_ERR_GROUP)DEVABORT("Null or invalid group");
    if(err!=MPI_SUCCESS)DEVABORT("Unclassified error");
#endif
    return currentCommunicator;
}

void MPIwrapper::print(string String){
    for(int n=0;n<Size();n++){
        if(Rank()==n)Out().MPIwrapper::Out::operator<<(String+"\n");
        Barrier();
    }
}

bool MPIwrapper::isMaster(MPI_Comm Comm){return Rank(Comm)==MPIwrapper_master;}
int  MPIwrapper::master(){return MPIwrapper_master;}
int  MPIwrapper::anySource(){
#ifndef _NOMPI_
    return MPI_ANY_SOURCE;
#else
    return 0;
#endif
}

MPIwrapper::File::File(const std::string Name)
{
    string name=Str(Name,"")+"_"+Rank();
    open(name.c_str(),std::ofstream::app);
}

int MPIwrapper::Rank(MPI_Comm Comm){
#ifndef _NOMPI_
    int r=0;
    MPI_Comm_rank(Comm,&r);
    return r;
#else
    return 0;
#endif
}

int MPIwrapper::Size(MPI_Comm Comm){
#ifndef _NOMPI_
    int s=1;
    MPI_Comm_size(Comm,&s);
    return s;
#else
    return 1;
#endif
}

void MPIwrapper::GatherV(std::complex<double> *Loc, int LSize, std::complex<double> *Glob, int *GSizes, int From){
#ifndef _NOMPI_
    MPISTART;
    if(Loc==Glob)ABORT("local and global storage must not be aliased");
    vector<int> disp(Size(),0);
    for(unsigned int k=1;k<disp.size();k++)disp[k]=disp[k-1]+GSizes[k-1];
    MPI_Gatherv(Loc,LSize,MPI_DOUBLE_COMPLEX,Glob,GSizes,disp.data(),MPI_DOUBLE_COMPLEX,From,MPIwrapper::communicator());
    MPISTOP;
#else
    std::memcpy(Glob,Loc,LSize*sizeof(*Glob));
#endif
}

void MPIwrapper::AllGatherV(std::complex<double> *Loc, int LSize, std::complex<double> *Glob, int *GSizes){
#ifndef _NOMPI_
    MPISTART;
    if(Loc==Glob)ABORT("local and global storage must not be aliased");
    vector<int> disp(Size(),0);
    for(unsigned int k=1;k<disp.size();k++)disp[k]=disp[k-1]+GSizes[k-1];
    MPI_Allgatherv(Loc,LSize,MPI_DOUBLE_COMPLEX,Glob,GSizes,disp.data(),MPI_DOUBLE_COMPLEX,MPIwrapper::communicator());
    MPISTOP;
#else
    memcpy(Glob,Loc,LSize*sizeof(*Glob));
#endif
}

void MPIwrapper::ScatterV(std::complex<double> *Glob, int *GSizes, std::complex<double> *Loc, int LSize, int From){
#ifndef _NOMPI_
    MPISTART;
    if(Loc==Glob)ABORT("local and global storage must not be aliased");
    vector<int> disp(Size(),0);
    for(unsigned int k=1;k<disp.size();k++)disp[k]=disp[k-1]+GSizes[k-1];
    // in case of aliasing, creat copy
    MPI_Scatterv(Glob,GSizes,disp.data(),MPI_DOUBLE_COMPLEX,Loc,LSize,MPI_DOUBLE_COMPLEX,From,MPIwrapper::communicator());
    MPISTOP;
#else
    memcpy(Loc,Glob,LSize*sizeof(*Glob));
#endif
}

void MPIwrapper::ISend(std::complex<double> *Data, int Size, int To, MPI_Request &Req,int Tag){
#ifndef _NOMPI_
    MPISTART;
    if(Tag==MPIwrapper::undefined)Tag=0;
    MPI_Isend(Data,Size,MPI_DOUBLE_COMPLEX,To,Tag,MPIwrapper::communicator(),&Req);
    MPISTOP;
#endif
}

void MPIwrapper::IRecv(std::complex<double> *Data, int Size, int From, MPI_Request &Req,int Tag){
#ifndef _NOMPI_
    MPISTART;
    if(Tag==MPIwrapper::undefined)Tag=0;
    MPI_Irecv(Data,Size,MPI_DOUBLE_COMPLEX,From,Tag,MPIwrapper::communicator(),&Req);
    MPISTOP;
#endif
}
void MPIwrapper::Recv(std::complex<double> *Data, int Size, int From, int Tag){
#ifndef _NOMPI_
    MPISTART;
    if(Tag==MPIwrapper::undefined)Tag=MPI_ANY_TAG;
    MPI_Recv(Data,Size,MPI_DOUBLE_COMPLEX,From,Tag,MPIwrapper::communicator(),MPI_STATUS_IGNORE);
    MPISTOP;
#endif
}
void MPIwrapper::Send(std::complex<double> *Data, int Size, int To, int Tag){
#ifndef _NOMPI_
    MPISTART;
    if(Tag==MPIwrapper::undefined)Tag=0;
    MPI_Send(Data,Size,MPI_DOUBLE_COMPLEX,To,Tag,MPIwrapper::communicator());
    MPISTOP;
#endif
}

void MPIwrapper::Barrier(MPI_Comm Comm){
#ifndef _NOMPI_
    MPISTART;
    MPI_Barrier(Comm);
    MPISTOP;
#endif
}

void MPIwrapper::Barrier_slowButHard(){
#ifndef _NOMPI_
    MPISTART;
    cout << flush;
    MPIwrapper::Barrier();
    int id = MPIwrapper::Rank();
    int nthreads = MPIwrapper::Size();
    int temp = 0;

    if(id != 0)
        MPIwrapper::Send(id,0,id);
    if(id == 0)
        for(int i=1;i<nthreads;i++)
            MPIwrapper::Recv(temp,i,i);

    if(id == 0)
        for(int i=1;i<nthreads;i++)
            MPIwrapper::Send(temp,i,i);
    if(id != 0)
        MPIwrapper::Recv(temp,0,id);
    cout << flush;
    MPIwrapper::Barrier();
    MPISTOP;
#endif
}


void MPIwrapper::Bcast(char* Data, int Size, int From){
#ifndef _NOMPI_
    MPISTART;
    MPI_Bcast(Data,Size,MPI_CHAR,From,MPIwrapper::communicator());
    MPISTOP;
#endif
}
void MPIwrapper::Bcast(int* Data, int Size, int From){
#ifndef _NOMPI_
    MPISTART;
    MPI_Bcast(Data,Size,MPI_INT,From,MPIwrapper::communicator());
    MPISTOP;
#endif
}

void MPIwrapper::Bcast(std::complex<double>* Data, int Size, int From){
#ifndef _NOMPI_
    MPISTART;
    MPI_Bcast(Data,Size,MPI_DOUBLE_COMPLEX,From,MPIwrapper::communicator());
    MPISTOP;
#endif
}

void MPIwrapper::Bcast(double* Data, int Size, int From){
#ifndef _NOMPI_
    MPISTART;
    MPI_Bcast(Data,Size,MPI_DOUBLE,From,MPIwrapper::communicator());
    MPISTOP;
#endif
}

void MPIwrapper::Bcast(string &Data, int From){
#ifndef _NOMPI_
    MPISTART;
    int cnt=Data.size();
    Bcast(&cnt,1,From);
    if(cnt>0){
        if(MPIwrapper::Rank()!=From)Data.assign(cnt,' ');
        Bcast(&Data[0],cnt,From);
    }
    else {
        Data.clear();
    }
    MPISTOP;
#endif
}

void MPIwrapper::Bcast(vector<string> &Data, int From){
#ifndef _NOMPI_
    MPISTART;
    //    if(Rank()==From){memcpy(msg,res.c_str(),res.size());}
    std::string all;
    std::vector<int> len;
    if(Rank()==From){
        for(unsigned int k=0;k<Data.size();k++){
            all+=Data[k];
            len.push_back(Data[k].length());
        }
        len.push_back(all.length());
    }
    // number of strings
    int cnt=len.size();
    Bcast(&cnt,1,From);

    // lengths of strings
    len.resize(cnt);
    Bcast(len.data(),len.size(),From);

    // all characters
    all.resize(len.back());
    Bcast(&all[0],len.back(),From);

    // reconstruct string vector
    Data.clear();
    len.pop_back();
    cnt=0;
    for(unsigned int k=0;k<len.size();k++){
        Data.push_back(all.substr(cnt,len[k]));
        cnt+=len[k];
    }
    MPISTOP;
#endif
}

void MPIwrapper::Send(int & x, int dest, int tag){
#ifndef _NOMPI_
    MPISTART;
    MPI_Send(&x,1,MPI_INT,dest,tag,MPIwrapper::communicator());
    MPISTOP;
#endif
}

void MPIwrapper::Send(int * x, int size, int dest, int Tag){
#ifndef _NOMPI_
    MPISTART;
    if(Tag==undefined)Tag=0;
    MPI_Send(x,size,MPI_INT,dest,Tag,MPIwrapper::communicator());
    MPISTOP;
#endif
}

void MPIwrapper::Send(char * x, int size, int dest, int Tag){
#ifndef _NOMPI_
    MPISTART;
    if(Tag==undefined)Tag=0;
    MPI_Send(x,size,MPI_CHAR,dest,Tag,MPIwrapper::communicator());
    MPISTOP;
#endif
}

void MPIwrapper::Send(std::string & Str, int dest){
#ifndef _NOMPI_
    MPISTART;
    char * c=const_cast<char*>(Str.data());
    int  siz=Str.size();
    MPIwrapper::Send(&siz,1,dest,0);
    MPIwrapper::Send(c,siz,dest,1);
    MPISTOP;
#endif
}

int MPIwrapper::Recv(std::string & Str, int dest){
#ifndef _NOMPI_
    MPISTART;
    int  siz;
    int from=MPIwrapper::Recv(&siz,1,dest,0);
    Str.resize(siz);
    char * c=const_cast<char*>(Str.data());
    MPIwrapper::Recv(c,siz,from,1);
    MPISTOP;
    return from;
#endif
    return 0;
}

void MPIwrapper::Recv(int & x, int dest, int tag){
#ifndef _NOMPI_
    MPISTART;
    // cout << "r: " << dest << " " << tag << endl;
    MPI_Recv(&x,1,MPI_INT,dest,tag,MPIwrapper::communicator(),MPI_STATUS_IGNORE);
    MPISTOP;
#endif
}

int MPIwrapper::Recv(int * x, int size, int dest, int Tag){
#ifndef _NOMPI_
    MPISTART;
    MPI_Status stat;
    if(Tag==undefined)Tag=MPI_ANY_TAG;
    MPI_Recv(x,size,MPI_INT,dest,Tag,MPIwrapper::communicator(),&stat);
    MPISTOP;
    return stat.MPI_SOURCE;
#endif
    return 0;
}

void MPIwrapper::Recv(char * x, int size, int dest, int Tag){
#ifndef _NOMPI_
    MPISTART;
    if(Tag==undefined)Tag=MPI_ANY_TAG;
    MPI_Recv(x,size,MPI_CHAR,dest,Tag,MPIwrapper::communicator(),MPI_STATUS_IGNORE);
    MPISTOP;
#endif
}


void MPIwrapper::Send(double & x, int dest, int tag){
#ifndef _NOMPI_
    MPISTART;
    MPI_Send(&x,1,MPI_DOUBLE,dest,tag,MPIwrapper::communicator());
    MPISTOP;
#endif
}

void MPIwrapper::Recv(double & x, int dest, int tag){
#ifndef _NOMPI_
    MPISTART;
    // cout << "r: " << dest << " " << tag << endl;
    MPI_Recv(&x,1,MPI_DOUBLE,dest,tag,MPIwrapper::communicator(),MPI_STATUS_IGNORE);
    MPISTOP;
#endif
}

//TIMER(allreduce1,)
//TIMER(allreduce2,)
//TIMER(allreduce3,)
//TIMER(allreduce4,)
//TIMER(allreduce5,)
//TIMER(allreduce6,)
//TIMER(allreduce7,)
void MPIwrapper::AllreduceMAX(double &send, double &recv){
#ifndef _NOMPI_
    MPISTART;
    //    STARTDEBUG(allreduce1);
    MPI_Allreduce(&send,&recv,1,MPI_DOUBLE,MPI_MAX,MPIwrapper::communicator());
    //    STOPDEBUG(allreduce1);
    MPISTOP;
#else
    recv = send;
#endif
}

void MPIwrapper::AllreduceSUM(complex<double> &send, complex<double> &recv){
#ifndef _NOMPI_
    MPISTART;
    //    STARTDEBUG(allreduce2);
    recv=0.;
    MPI_Allreduce(&send,&recv,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPIwrapper::communicator());
    //    STOPDEBUG(allreduce2);
    MPISTOP;
#else
    recv = send;
#endif
}

void MPIwrapper::AllreduceSUM(complex<double> * Data, int Siz){
#ifndef _NOMPI_
    if(Size()==1)return;
    MPISTART;
    vector<complex<double> > dat(Siz,0.);
    MPI_Allreduce(Data,dat.data(),Siz,MPI_DOUBLE_COMPLEX,MPI_SUM,MPIwrapper::communicator());
    for(unsigned int k=0;k<dat.size();k++)Data[k]=dat[k];
    MPISTOP;
#endif
}

void MPIwrapper::AllreduceSUM(int * Data, int Siz){
#ifndef _NOMPI_
    MPISTART;
    vector<int> dat(Siz);
    MPI_Allreduce(Data,dat.data(),Siz,MPI_INT,MPI_SUM,MPIwrapper::communicator());
    for(unsigned int k=0;k<dat.size();k++)Data[k]=dat[k];
    MPISTOP;
#endif
}

void MPIwrapper::AllreduceSUM(double * Data, int Siz){
#ifndef _NOMPI_
    MPISTART;
    vector<double> dat(Siz);
    MPI_Allreduce(Data,dat.data(),Siz,MPI_DOUBLE,MPI_SUM,MPIwrapper::communicator());
    for(unsigned int k=0;k<dat.size();k++)Data[k]=dat[k];
    MPISTOP;
#endif
}

void MPIwrapper::AllreduceMAX(int * Data, int Siz){
#ifndef _NOMPI_
    MPISTART;
    vector<int> dat(Siz);
    MPI_Allreduce(Data,dat.data(),Siz,MPI_INT,MPI_MAX,MPIwrapper::communicator());
    for(unsigned int k=0;k<dat.size();k++)Data[k]=dat[k];
    MPISTOP;
#endif
}

void MPIwrapper::AllreduceMAX(double * Data, int Siz){
#ifndef _NOMPI_
    MPISTART;
    vector<double> dat(Siz);
    MPI_Allreduce(Data,dat.data(),Siz,MPI_DOUBLE,MPI_MAX,MPIwrapper::communicator());
    for(unsigned int k=0;k<dat.size();k++)Data[k]=dat[k];
    MPISTOP;
#endif
}

void MPIwrapper::AllreduceMAX(long long &Data){
#ifndef _NOMPI_
    MPISTART;
    long long data;
    MPI_Allreduce(&Data,&data,1,MPI_LONG_LONG_INT,MPI_MAX,MPIwrapper::communicator());
    data=Data;
    MPISTOP;
#endif
}

string MPIwrapper::Get_processor_name(){
#ifdef _NOMPI_
    return "0";
#else
    MPISTART;
    char* name = new char[MPI_MAX_PROCESSOR_NAME];
    int length = 0;
    MPI_Get_processor_name(name,&length);
    string aNiceString(name, length);
    delete[] name;
    return aNiceString;
    MPISTOP;
#endif
}

static bool _isMPI=false;
bool MPIwrapper::isMPI(){return _isMPI;}

static Str strRank;
void MPIwrapper::Init(int argc, char *argv[]) {
#ifndef _NOMPI_
    MPISTART;
    _isMPI=true;
    char** mpiArgv=const_cast<char**>(argv);
    MPI_Init(&argc,&mpiArgv);
    MPISTOP;
    strRank="<"+tools::str(Rank())+">";
#endif
}

void MPIwrapper::Finalize() {
#ifndef _NOMPI_
    _isMPI=false;
    MPI_Finalize();
#endif
}

void MPIwrapper::Exit(int i) {
    if(Rank()==0) cout << "exiting with code " << i << endl;
#ifndef _NOMPI_
    Barrier();
    Finalize();
    sleep(1);
#endif
    if(i==0)  exit(i);
    else abort();
}

void MPIwrapper::printAll(string S){
    for (int k=0;k<Size();k++){
        string s=Str("<","")+Rank()+">"+S;
        Bcast(s,k);
        if(isMaster())std::cout<<s<<std::endl;
    }
}

void MPIwrapper::Abort(int i){
#ifndef _NOMPI_
    sleep(1);
#endif
    if(i==0){
        Finalize();
        exit(0);
    }
    else
#ifndef _NOMPI_
        MPI_Abort(MPI_COMM_WORLD,i);
#else
        abort();
#endif
}

void MPIwrapper::Wait(MPI_Request *r){
#ifndef _NOMPI_
    MPISTART;
    MPI_Wait(r,MPI_STATUS_IGNORE);
    MPISTOP;
#endif
}

void MPIwrapper::Waitall(vector<MPI_Request> &array_of_requests){
#ifndef _NOMPI_
    if(array_of_requests.size() != 0){
        MPISTART;
        MPI_Waitall(array_of_requests.size(),array_of_requests.data(),MPI_STATUS_IGNORE);
        MPISTOP;
    }
#endif
}
void MPIwrapper::Waitall(const vector<Buffer> & Buf){
#ifndef _NOMPI_
    MPISTART;
    vector<MPI_Request> req;
    for(unsigned int k=0;k<Buf.size();k++)
        if(Buf[k].size()>0)req.push_back(Buf[k].req);
    Waitall(req);
    MPISTOP;
#endif
}

void MPIwrapper::ISend(vector<complex<double> > &coef, int dest, int tag, MPI_Request *r){
#ifndef _NOMPI_
    MPISTART;
    MPI_Isend(coef.data(),coef.size(),MPI_DOUBLE_COMPLEX,dest,tag,MPIwrapper::communicator(),r);
    MPISTOP;
#endif
}

void MPIwrapper::IRecv(vector<complex<double> > &coef, int dest, int tag, MPI_Request *r){
#ifndef _NOMPI_
    MPISTART;
    MPI_Irecv(coef.data(),coef.size(),MPI_DOUBLE_COMPLEX,dest,tag,MPIwrapper::communicator(),r);
    MPISTOP;
#endif
}

void MPIwrapper::ISend(int &x, int dest, int tag, MPI_Request *r){
#ifndef _NOMPI_
    MPISTART;
    MPI_Isend(&x,1,MPI_INT,dest,tag,MPIwrapper::communicator(),r);
    MPISTOP;
#endif
}

void MPIwrapper::IRecv(int &x, int dest, int tag, MPI_Request *r){
#ifndef _NOMPI_
    MPISTART;
    MPI_Irecv(&x,1,MPI_INT,dest,tag,MPIwrapper::communicator(),r);
    MPISTOP;
#endif
}

void MPIwrapper::ISend(vector<int> &x, int dest, int tag, MPI_Request *r){
#ifndef _NOMPI_
    MPISTART;
    MPI_Isend(x.data(),x.size(),MPI_INT,dest,tag,MPIwrapper::communicator(),r);
    MPISTOP;
#endif
}

void MPIwrapper::IRecv(vector<int> &x, int dest, int tag, MPI_Request *r){
#ifndef _NOMPI_
    MPISTART;
    MPI_Irecv(x.data(),x.size(),MPI_INT,dest,tag,MPIwrapper::communicator(),r);
    MPISTOP;
#endif
}
