// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "wavefunction.h"

//resolve forward declarations
#include "threads.h"
#include "coefficients.h"
#include "discretization.h"
#include "discretizationDerived.h"
//#include "operator.h"
#include "discretizationGrid.h"
#include "index.h"

using namespace std;
using namespace tools;

Wavefunction::Wavefunction(const Discretization* D, double Time) :time(Time),_map(false){coefs=new Coefficients(D->idx());}
Wavefunction::Wavefunction(double Time, Coefficients *C):time(Time),_map(false){coefs=new Coefficients(*C);}
Wavefunction::Wavefunction(const Index* Idx, double Time):time(Time),_map(false){coefs=new Coefficients(Idx);}

Wavefunction::~Wavefunction(){if(not _map){delete coefs; coefs=0;}}


Wavefunction& Wavefunction::operator=(const Wavefunction& rhs)
{
    if(&rhs==this)return *this;

    if(coefs==0){
        coefs=new Coefficients(*rhs.coefs);
        coefs->treeOrderStorage();
    }
    else{
        (*coefs)=(*(rhs.coefs));
    }
    time = rhs.time;
    return *this;
}

Wavefunction::Wavefunction(const Wavefunction& rhs)
{
    coefs=new Coefficients(*rhs.coefs);
    coefs->treeOrderStorage();
    time=rhs.time;
}

const Wavefunction Wavefunction::map(double Time, const Coefficients* C)
{
    Wavefunction wf;
    (&wf)->time=Time;
    (&wf)->coefs=const_cast<Coefficients*>(C);
    return wf;
}

Wavefunction& Wavefunction::operator+=(const Wavefunction& rhs)
{
    (*coefs)+=(*(rhs.coefs));
    return *this;
}

Wavefunction& Wavefunction::operator-=(const Wavefunction& rhs)
{
    (*coefs) -= (*(rhs.coefs));
    return *this;
}

Wavefunction& Wavefunction::operator*=(std::complex< double > x)
{
    (*coefs) *= x;
    return *this;
}

void Wavefunction::setToZero() {
    coefs->setToZero();
}

void Wavefunction::setToConstant(complex<double> c){
    coefs->setToConstant(c);
}


void Wavefunction::setToRandom(int seed){
    if(seed!=0)srand48(seed);
    coefs->setToRandom();
}

void Wavefunction::makeContinuous() {
    coefs->makeContinuous();
}
void Wavefunction::show(int Digits){cout<<coefs->str(Digits)<<endl;}

/// header: reposition file to beginning and write header information
void Wavefunction::write(ofstream &stream, bool header) const{
    coefs->write(stream,header);
    tools::write(stream,time);
    stream.flush(); // this is important, as else the outfile may remain dangling
}
/// header: reposition file to beginning and write header information
void Wavefunction::print(ofstream &stream, string Header) const {
    if(MPIwrapper::isMaster(MPIwrapper::communicator())){
        if(Header!=""){
            stream <<"# "+Header<<endl;
            stream <<"# time abs(0), arg(0), ..."<<endl;
        }
        stream<<time;
    }
    coefs->print(stream);
}
/// header: reposition file to beginning and read header information
void Wavefunction::read(ifstream &stream,bool header){
    ///< header: read from beginning of file, including structural information
    coefs->read(stream,header);
    tools::read(stream,time);
}

double Wavefunction::readLastTime(ifstream &Inp){
    ABORT("havn't gotten this to work...");
    if(not Inp.is_open())ABORT("not a good input stream");
    double t=-DBL_MAX;
    Inp.seekg (8,Inp.end);
    cout<<" "<<t<<" "<<Inp.tellg()<<endl;

    tools::read(Inp,t);
    return t;
}

double Wavefunction::maxCoeff()
{
    return coefs->maxCoeff();
}
