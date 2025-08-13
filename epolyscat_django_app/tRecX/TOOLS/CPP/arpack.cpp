// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "mpiWrapper.h"
#include "arpack.h"
#include "tools.h"
#include "constants.h"
#include "printOutput.h"

#include <vector>

using namespace std;

// rather brute force definition of the interfaces to the arpack F77 routines (may not be very portable)
extern "C"{

void znaupd_(int*,char*,int*,char*,int*,double*,complex<double>*,int*,complex<double>*,int*,
             int*,int*,complex<double>*,complex<double>*,int*,double*,int*);

void zneupd_(int*,char*,int*,complex<double>*,complex<double>*,int*,complex<double>*,complex<double>*,
             char*,int*,char*,int*,double*,complex<double>*,int*,complex<double>*,int*,
             int*,int*,complex<double>*,complex<double>*,int*,double*,int*);
}

map<string,string> Arpack::InterfaceF77::whichList;
void Arpack::InterfaceF77::whichListSet(){
    whichList["SmallReal"]="SR";
    whichList["LargeReal"]="LR";
    whichList["SmallImag"]="SI";
    whichList["LargeImag"]="LI";
    whichList["LargeAbs"]="LM";
    whichList["SmallAbs"]="SM";
}

void Arpack::eigenValues(vector<complex<double> > & Eval,unsigned int Nvec,const string & Which)
{
    _I_F77.iter(this,Nvec,Which,false,false);
    Eval=eval;
}

void Arpack::eigenStdVec(vector<complex<double> > & Eval,vector<vector<complex<double> > >& Rvec,unsigned int Nvec,const string & Which,bool Restart)
{
    _I_F77.iter(this,Nvec,Which,Restart,true);
    Rvec.clear();
    Eval=eval;
    for(unsigned int k=0;k<eval.size();k++){
        Rvec.push_back(std::vector<std::complex<double> >(0));
        for (unsigned int i=0;i<_lvec;i++)Rvec[k].push_back(rvec[i+k*_lvec]);
    }
}

void Arpack::verify() {
    vector<complex<double> > zero(_lvec);
    double maxErr=0.,minErr=DBL_MAX;
    for(unsigned int k=0;k<eval.size();k++){
        apply(rvec.data()+k*_lvec,zero.data());
        double norm=0.;
        for (unsigned int i=0;i<_lvec;i++){
            norm=max(norm,abs(zero[i])); // L1 norm
            zero[i]-=eval[k]*rvec[i+k*_lvec];
        }
        double err=0.;
        for(unsigned int i=0;i<_lvec;i++)err=max(err,abs(zero[k])/norm); // L1 error
        minErr=min(minErr,err);
        maxErr=max(maxErr,err);
    }
    PrintOutput::message("verified "+tools::str(int(eval.size()))+" eigenvectors of length "+tools::str(_lvec)
                         +", min/max error ="+tools::str(minErr,3)+"/"+tools::str(maxErr,3));
}

void Arpack::eigenIter(unsigned int Nvec, std::string Which,bool Restart,bool ComputeVectors){
    _I_F77.iter(this,Nvec,Which,Restart,ComputeVectors);
}

void Arpack::InterfaceF77::iter(Arpack* Arp,int Nvec, std::string Which, bool Restart, bool computeVectors, bool Parallel){

    if(Nvec+1>int(Arp->_lvec))ABORT(Sstr+"at most N-1 eigenvalues can be calculated, got"+Nvec);

    // adjust storage
    Arp->eval.resize(Nvec+1); // used for temporary storage by zneupd, will be truncted to size=Nvec
    if(computeVectors)Arp->rvec.resize(Nvec*Arp->_lvec);

    if(whichList[Which]=="")ABORT("illegal Which="+Which+", admissible values: "+tools::listMapKeys(whichList,",",&Which));

    // initialze the communication parameters
    vector<int> iparam(11);
    iparam[0]=1;       // exact shifts
    iparam[2]=_maxIter;// maximum number of iterations
    iparam[3]=1;       // must be =1
    iparam[6]=1;       // solve standard problem

    // decide restart, put starting vector if needed
    int info;
    vector<complex<double> >resid(Arp->_lvec);
    if(Restart){
        if(Arp->rvec.size()<Arp->_lvec)ABORT("must supply initial vector for Restart");
        for(unsigned int i=0;i<Arp->_lvec;i++)resid[i]=Arp->rvec[i];
        info=1;
        PrintOutput::message("Arpack used with restart option, criterion = "+Which,0);
    }
    else
        info=0;

    int vecs=0;
    if(computeVectors)vecs=1;

    // auxiliary variables for znaupd
    int ido=0,lvec=Arp->_lvec,nvec=Nvec,ncv=min(int(2*Nvec+5),int(Arp->_lvec-1));
    char bmat='I';
    // convert string to char*
    vector<char> which(whichList[Which].begin(),whichList[Which].end());
    which.push_back('\0');
    vector<complex<double> > v, workd, workl;
    vector<double> rwork;

    if(_seriel or MPIwrapper::isMaster()){
        v.resize(Arp->_lvec*ncv);
        workl.resize(3*ncv*ncv+5*ncv);
        rwork.resize(ncv);
    }
    workd.resize(3*Arp->_lvec);

    vector<int> ipntr(14);
    int lworkl=workl.size();

    while (ido!=99) {
        if(_seriel or MPIwrapper::isMaster()){
            znaupd_(&ido,&bmat,&lvec,which.data(),&nvec,&_tolerance,resid.data(),&ncv,v.data(),&lvec,
                    iparam.data(),ipntr.data(),workd.data(),workl.data(),&lworkl,rwork.data(),&info);
            if(info!=0){
                if(info!=1)ABORT("znaupd failed, info="+tools::str(info));
                ido=99;// info=1: converged to maxIter
            }
        }

        if(not _seriel)MPIwrapper::Bcast(&ido,1,MPIwrapper::master());
        if(not _seriel)MPIwrapper::Bcast(ipntr.data(),ipntr.size(),MPIwrapper::master());
        switch (ido){
        default: ABORT("undefined value of ido: "+tools::str(ido));
        case 1:
        case -1:
            Arp->apply(workd.data()+ipntr[0]-1,workd.data()+ipntr[1]-1);
            break;
        case 99:
            // additional auxiliary variables for zneupd
            char howmny='A';
            vector<int> select_v(ncv);
            complex<double>sigma;
            vector<complex<double> > workev(2*ncv);
            if(_seriel or MPIwrapper::isMaster()){
                zneupd_(&vecs,&howmny,select_v.data(),Arp->eval.data(),Arp->rvec.data(),&lvec,&sigma,
                        workev.data(),&bmat,&lvec,which.data(),&nvec,&_tolerance,resid.data(),
                        &ncv,v.data(),&lvec,iparam.data(),ipntr.data(),workd.data(),
                        workl.data(),&lworkl,rwork.data(),&info);
                if(info!=0)ABORT("zneupd failed, info="+tools::str(info));
            }
            Arp->eval.pop_back(); // original size is Nvec+1
            if(not _seriel)MPIwrapper::Bcast(Arp->eval.data(),Arp->eval.size(),MPIwrapper::master());
        }
    }
}

void ArpackStandard::apply(const complex<double> *X, complex<double> *Y) {
    for(complex<double> *y=Y;y<Y+_lvec;y++)*y=0;
    unsigned int m=0;
    for(unsigned int i=0;i<_lvec;i++)
        for(complex<double> *y=Y;y<Y+_lvec;y++,m++)
            *y+=mat[m]*(*(X+i));
}

void ArpackStandard::test(){
    // generate a random matrix
    vector<complex<double> > mat,eval;
    vector<vector<complex<double> > >  evec;
    const unsigned int N=500;

    /* initialize random seed: */
    srand (1);
    double norm=2.*math::pi/double(RAND_MAX);
    for(unsigned int k=0;k<N*N;k++)
        mat.push_back(complex<double>(rand()*norm,rand()*norm));

    // create the ArpackStandard
    ArpackStandard A(mat,1.e-10,100);
    A.eigenStdVec(eval,evec,5,"SmallReal");
    A.verify();
}
