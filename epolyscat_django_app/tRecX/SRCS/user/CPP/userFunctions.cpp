// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "userFunctions.h"
#include "tools.h"

using namespace std;


/// all functions must be added to the function table
void userFunctionsList(){
    basisMatFuncSet(new basisMatUserRedefined());
    basisMatFuncSet(new basisMatPyrA(00));
    basisMatFuncSet(new basisMatPyrA(01));
    basisMatFuncSet(new basisMatPyrA(10));
    basisMatFuncSet(new basisMatPyrA(11));
    basisMatFuncSet(new basisMatPyrA(33));
    basisMatFuncSet(new basisMatPyrA(44));
    basisMatFuncSet(new basisMatHarMorse(0));
    basisMatFuncSet(new basisMatHarMorse(1));
    basisMatFuncSet(new basisMat2Dhelium());
    basisMatFuncSet(new basisMatXtimesY());
}

basisMatUserRedefined::basisMatUserRedefined():basisMatFunc("user"){}

static bool firstUserRedefined=true;

std::complex<double> basisMatUserRedefined::operator()(std::vector<std::complex<double> > X) const
{
    std::complex<double> val;

    if(firstUserRedefined){
        firstUserRedefined=false;
        PrintOutput::warning(
                "user{...} does not have fixed meaning - check code in "
                +string(__FILE__)+" at line "+tools::str(__LINE__+3));
    }
    //====== change code between these lines only ====================

    // check number of arguments for functions
    basisMatFunc::needSize(X,2);

    // evaluate function
    val=0.5*(pow(X[0],2)+pow(X[1],2));

    //================================================================

    return val;
}



double ev2hartree=1./27.21138505;

double rmn=2.78;
double mm=4.*21894.16+4.*1837.36;
double mn=25532.65;
double mh=1837.36;


double De1=4.979*ev2hartree;
double r1=1.927;
double a1=1.137;
double wc1=0.1096*ev2hartree;

double E02=4.805*ev2hartree;
double De21=4.979*ev2hartree;
double r21=1.882;
double a21=1.293;

double l22=1.248*ev2hartree;

double Ab22=2.644*ev2hartree;
double De22=3.956*ev2hartree;
double r22=2.216;
double a22=1.325;
double wc2=0.1096*ev2hartree;

double l12max=0.237*ev2hartree;
double d12= 3.679;
double beta12=1.369;

double x=rmn*(mm/(mn+mm));
complex<double> rnh(complex<double> r, double Q){
    return r-x;
}

complex<double> l12(complex<double> r, double Q){
    return 0.5*l12max*(1.-tanh((rnh(r,Q)-d12)/beta12));
}

complex<double> nu11(complex<double> r, double Q){
    return De1*pow(1.-exp(-a1*(rnh(r,Q)-r1)),2);
}

complex<double> nu21(complex<double> r, double Q){
    return De21*pow(1.-exp(-a21*(rnh(r,Q)-r21)),2)+E02;
}

complex<double> nu22(complex<double> r, double Q){
    return Ab22*exp(-a22*(rnh(r,Q)-r22))+De22;
}

complex<double> V11(complex<double> r,double Q){
    return nu11(r,Q)+0.5*wc1*pow(Q,2);
}

complex<double> V22(complex<double> r,double Q){
    return 0.5*(nu21(r,Q)+nu22(r,Q))-0.5*sqrt(pow(nu21(r,Q)-nu22(r,Q),2)+4.*pow(l22,2))+0.5*wc2*pow(Q,2);
}

complex<double> V12(complex<double> r,double Q){
    return l12(r,Q)*Q;
}

complex<double> V1a(complex<double> r,double Q){
    return 0.5*(V11(r,Q)+V22(r,Q))-0.5*sqrt(pow(V11(r,Q)-V22(r,Q),2)+4.*pow(l12(r,Q)*Q,2));
}

complex<double> V2a(complex<double> r,double Q){
    return 0.5*(V11(r,Q)+V22(r,Q))+0.5*sqrt(pow(V11(r,Q)-V22(r,Q),2)+4.*pow(l12(r,Q)*Q,2));
}

basisMatPyrA::basisMatPyrA(int Kind):basisMatFunc("PyrA["+tools::str(Kind,2,'0')+"]"),kind(Kind){}
std::complex<double> basisMatPyrA::operator()(std::vector<std::complex<double> >X) const
{
    basisMatFunc::needSize(X,2); // check whether called with correct number of arguments
    switch (kind){
    case 00: return V11(X[1],X[0].real());
    case 01: return V12(X[1],X[0].real());
    case 10: return V12(X[1],X[0].real());
    case 11: return V22(X[1],X[0].real());
    case 33: return V1a(X[1],X[0].real());
    case 44: return V2a(X[1],X[0].real());
    default: ABORT("undefined kind="+tools::str(kind));
    }
}
basisMatHarMorse::~basisMatHarMorse(){delete morse0; delete morse1;}
basisMatHarMorse::basisMatHarMorse(int Kind):basisMatFunc("harMorse["+tools::str(Kind,1)+"]"),kind(Kind)
{
    morse0= new AlgebraMorse("Morse[0.5,2]");
    morse1= new AlgebraMorse("Morse[0.25,4]");
}
std::complex<double> basisMatHarMorse::operator()(std::vector<std::complex<double> >X) const
{
    basisMatFunc::needSize(X,2); // check whether called with correct number of arguments
    switch (kind){
    case 0:
        // <0.5*pow[2](Q)-1><1>+<1><Morse[0.5,2]>)
        return 0.5*X[0]*X[0]-1.+morse0->val(X[1]);
    case 1:
        //+<0.125.pow[2](Q)><1>+0.01<1><Morse[1.,4]>
        // <0.125*pow[2](Q)><1>+0.01<1><Morse[0.25,4]>
        return 0.125*X[0]*X[0]+0.01*morse1->val(X[1]);

    default: ABORT("undefined kind="+tools::str(kind));
    }
}

basisMatXtimesY::basisMatXtimesY():basisMatFunc("XtimesY"){}

std::complex<double> basisMatXtimesY::operator()(std::vector<std::complex<double> > X) const
{
    basisMatFunc::needSize(X,2); // check whether called with correct number of arguments
    return X[0]*X[1];
}


basisMat2Dhelium::basisMat2Dhelium():basisMatFunc("EERepulsion2D"){}

std::complex<double> basisMat2Dhelium::operator()(std::vector<std::complex<double> > X) const
{
    basisMatFunc::needSize(X,2); // check whether called with correct number of arguments
    return 1./std::sqrt(std::pow(X[0]-X[1],2)+0.3);
}

