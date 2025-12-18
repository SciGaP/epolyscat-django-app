// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "schroedingerMachine.h"

#include <string>
#include <complex>
#include "readInput.h"
#include "printOutput.h"

#include "operatorTree.h"
#include "operatorDefinition.h"
#include "basisDvr.h"
#include "index.h"
#include "coefficients.h"
#include "eigenSolver.h"
#include "str.h"
#include "polylagrange.h"

std::string SchroedingerMachine::_task;
std::string SchroedingerMachine::_potential;
bool SchroedingerMachine::_exactIntegrals=true;

static std::vector<const BasisDVR*> _trainBas;

static SchroedingerMachine* _main;
static SchroedingerMachine& m(){return * _main;}

static Eigen::MatrixXd _weigBasic; //stub
static Eigen::MatrixXd _nodeBasic; //stub
void setBasic(int Nelem1){
    // at present, creates Nelem-1 identical elements
    std::vector<double> node0,weig0;
    // lobatto quadrature on [0,1]
    OrthogonalLegendre().quadratureLobatto(_trainBas[0]->order(),node0,weig0);
    for(double &x: node0)x=0.5*(x+1.);
    for(double &x: weig0)x=0.5*x;

    _nodeBasic.resize(node0.size(),Nelem1-1);
    _weigBasic.resize(node0.size(),Nelem1-1);
    for(int n=0;n<_weigBasic.cols();n++){
        _nodeBasic.col(n)=Eigen::Map<Eigen::VectorXd>(node0.data(),node0.size());
        _weigBasic.col(n)=Eigen::Map<Eigen::VectorXd>(weig0.data(),weig0.size());
    }
}

static bool intervalsSharePoint(double A0,double A1,double B0,double B1,double Epsilon){
    return A1>B0-Epsilon and B1>A0-Epsilon;
}

SchroedingerMachine::SchroedingerMachine(const Index *Idx, std::string Potential)
{
    _ham.reset(new OperatorTree("Ham",OperatorDefinition("0.5<<Laplacian>>+"+Potential,Idx->hierarchy()),Idx,Idx));

    EigenSolver es(DBL_MAX,DBL_MAX,false);
    es.select("SmallReal[1]");
    es.compute(_ham.get());
    PrintOutput::message(Sstr+"exact ground state energy: "+es.eigenvalues()[0]);
}

void SchroedingerMachine::read(){
    std::string task;
    ReadInput::main.read("SchroedingerMachine","task",_task,"NONE","machine learning, options: ground state...find ground state for the given potential");
    ReadInput::main.read("SchroedingerMachine","potential",_potential,"<<Harmonic>>","machine learning, options: ground state...find ground state for the given potential");
}


bool SchroedingerMachine::exactIntegrals(bool ExactIntegrals){
    if(_task!="NONE")return _exactIntegrals; return ExactIntegrals;
}

int sm_NElemTrain(){return _trainBas.size();}

void sm_SizesTrain(const char *File, int *NPsi, int *SizePsi){

    std::ifstream inStr(File);
    // get length of file:
    inStr.seekg (0, inStr.end);
    int lenIn=inStr.tellg();
    Index* idx=new Index(inStr,true); // rewind and read Index

    Coefficients c(idx);
    _trainBas.clear();

    for(const Coefficients* f=c.firstFloor();f!=0;f=f->nodeRight())
        _trainBas.push_back(dynamic_cast<const BasisDVR*>(f->idx()->basis()));
    for(auto b: _trainBas)
        if(b->order()!=_trainBas.front()->order())
            DEVABORT("for now, only constant element order allowed");

    *SizePsi=2*_trainBas.front()->order()*_trainBas.size();
    lenIn-=inStr.tellg();
    *NPsi=lenIn/((idx->size()*2+1)*sizeof(double));
    delete idx;
}

void sm_Quadrature(int Order, int Nelem, double *Node0, double *Weig0){
    setBasic(Nelem+1);
    if(_nodeBasic.rows()!=Order)DEVABORT("quadrature orders do not match");
    for(int in=0;in<Order*Nelem;in++){
        Node0[in]=_nodeBasic.data()[in];
        Weig0[in]=_weigBasic.data()[in];
    }
}

void sm_FuncPsiC(int Order, int Nelem, const double *X, const double *PsiIn, double *PsiOut, const double *Grad){
    SchroedingerMachine sm;
    if(Grad)
        Eigen::Map<Eigen::MatrixXd>(PsiOut,Order,Nelem)=
                sm.funcPhi(Eigen::Map<const Eigen::MatrixXd>(X,Order,Nelem),
                           Eigen::Map<const Eigen::MatrixXcd>(reinterpret_cast<const std::complex<double>*>(PsiIn),Order,20),
                           Eigen::Map<const Eigen::MatrixXcd>(reinterpret_cast<const std::complex<double>*>(Grad),Order,Nelem));
    else
        Eigen::Map<Eigen::MatrixXd>(PsiOut,Order*2,Nelem)=
                sm.funcPhi(Eigen::Map<const Eigen::MatrixXd>(X,Order,Nelem),
                           Eigen::Map<const Eigen::MatrixXcd>(reinterpret_cast<const std::complex<double>*>(PsiIn),Order,20),
                           Eigen::MatrixXcd());
}




void sm_LossL2(int Nelem, const double *Elem, const double *Psi, double* Res, const double* Grad){
    SchroedingerMachine sm;
    Eigen::VectorXf elem(Nelem);
    Eigen::VectorXf grad;
    for(int k=0;k<Nelem;k++)elem(k)=Elem[k];
    if(Grad){
        grad.resize(Nelem);
        for(int k=0;k<Nelem;k++)grad(k)=Grad[k];
    }
    Eigen::Map<Eigen::MatrixXd>(Res,1,1)=sm.LossL2(elem,Eigen::Map<const Eigen::MatrixXcd>(reinterpret_cast<const std::complex<double>*>(Psi),10,20),grad);
}

static PolyLagrange<double> _bas;
static std::vector<double> _basWeig;
void sm_DataTrain(const char *File, double *TrainData){
    std::cout<<"inStr "<<File<<std::endl;
    std::ifstream inStr(File);
    Index* idx=new Index(inStr,true); // re-read index
    Coefficients c(idx);

    double time(0.);
    double * dat=TrainData;

    while(c.read(inStr,false)){
        if(_trainBas.front()->nBeg()==1){
            // left dirchlet boundary condition
            dat[0]=0;
            dat[1]=0;
            dat+=2;
        }
        tools::read(inStr,time);
        for(std::complex<double> *d=c.data();d<c.data()+c.size();d++,dat+=2){
            dat[0]=d->real();
            dat[1]=d->imag();
        }
        if(_trainBas.back()->size()!=_trainBas.back()->order()){
            // right dirchlet boundary condition
            dat[0]=0;
            dat[1]=0;
            dat+=2;
        }
    }
    delete idx;

    // basic lagrange polynomial
    std::vector<double> node,weig;
    OrthogonalLegendre().quadratureLobatto(_trainBas.front()->order(),node,weig);
    for(auto &q: node)q=0.5*q+0.5;
    for(auto &q: weig)q=0.5*q;
    _bas=PolyLagrange<double>(node);
    _basWeig=weig;



    Sout+"OK loaded"+File+TrainData[0]+TrainData[1]+Sendl;
}

void sm_StrucTrain(int *SizeElem, double *BElem){
    for(int n=0;n<_trainBas.size();n++){
        SizeElem[n]=_trainBas[n]->size();
        BElem[n]=_trainBas[n]->lowBound();
    }
    BElem[_trainBas.size()]=_trainBas.back()->upBound();
}


Eigen::MatrixXd & weigBasic(){return _weigBasic;}
Eigen::MatrixXd & nodeBasic(){return _nodeBasic;}


Eigen::VectorXd SchroedingerMachine::funcC(const Eigen::VectorXf &Elem,const Eigen::VectorXf Grad) const{
    Eigen::VectorXd c=Eigen::VectorXd::Zero(Elem.size()-1);
    if(Grad.size()){
        c(0)=Grad(0)/Elem(0);
        for(int k=1;k<c.size();k++)c(k)=c(k-1)+Grad(k+1);
    }
    else {
        c(0)=log(Elem(0));
        for(int k=1;k<c.size();k++)c(k)=c(k-1)+Elem(k+1);
    }
    return c;
}

Eigen::MatrixXd SchroedingerMachine::funcX(const Eigen::VectorXf &Elem,const Eigen::VectorXf Grad) const{
    Eigen::MatrixXd x(_nodeBasic);
    Eigen::VectorXd c(funcC(Elem,Grad));
    for(int n=0;n<x.cols();n++)
        Grad.size()?x.col(n)=x.col(n)*Grad(n+1)+Eigen::VectorXd::Constant(x.rows(),c(n))
                :   x.col(n)=x.col(n)*Elem(n+1)+Eigen::VectorXd::Constant(x.rows(),c(n));
    return x;
}

Eigen::MatrixXd SchroedingerMachine::funcW(const Eigen::VectorXf &Elem,const Eigen::VectorXf Grad) const{
    Eigen::MatrixXd w(_weigBasic);
    if(Grad.size())for(int n=0;n<w.cols();n++)w.col(n)*=Grad(n+1);
    else           for(int n=0;n<w.cols();n++)w.col(n)*=Elem(n+1); // at present, unfortunate numbering of Elem
    return w;
}

static Eigen::MatrixXi ppos;
Eigen::MatrixXd SchroedingerMachine::funcPhi(const Eigen::MatrixXd &X, const Eigen::MatrixXcd &Psi, const Eigen::MatrixXcd Grad) const{

    // here the core of the algorithm
    Eigen::MatrixXcd psiX=Eigen::MatrixXcd::Zero(X.rows(),X.cols());
    Eigen::MatrixXi epos=Eigen::MatrixXi::Constant(X.rows(),X.cols(),-1);
    int m=0;
    std::vector<std::complex<double>>v,w;
    // assumes X and _trainBas are sorted
    for(int n=0;n<X.cols();n++){
        for(int i=0;i<X.rows();i++){
            // locate X in training element
            if(X(i,n)<_trainBas[m]->lowBound())continue;
            if(X(i,n)>_trainBas[m]->upBound()){
                if(X(i,n)>_trainBas.back()->upBound())break;
                while(X(i,n)>_trainBas[m]->upBound())m++;
                v.clear();
            }

            double x=X(i,n);
            if(abs(x-_trainBas[m]->lowBound())<1.e-3)x=_trainBas[m]->lowBound();
            if(abs(x-_trainBas[m]->upBound() )<1.e-3)x=_trainBas[m]->upBound();

            if(Grad.size())v=(_trainBas[m]->der(x));
            else           v=(_trainBas[m]->val(x));

            // supplement first or last by Dirichlet zeros
            if(m==0)v.insert(v.begin(),0.);
            if(m+1==_trainBas.size())v.push_back(0.);
            psiX(i,n)=Eigen::Map<Eigen::VectorXcd>(v.data(),v.size()).dot(Psi.col(m));
            if(epos(i,n)==-1)epos(i,n)=m;
            else epos(i,n)=-m;

        }
    }
    Eigen::MatrixXd res=Eigen::Map<Eigen::MatrixXd>(reinterpret_cast<double*>(psiX.data()),2*X.rows(),X.cols());
//    if(ppos.size()==0)ppos=epos;
//    if(not (ppos-epos).isZero()){
//        std::cout<<"epos\n"<<epos<<std::endl;
//        std::cout<<"diff\n"<<(epos-ppos)<<std::endl;
//        std::cout<<"x "<<X(0,0)<<" "<<X(X.rows()-1,X.cols()-1)<<std::endl;
//        if(epos(1,0)==-1){
//           std::cout<<"X\n"<<X<<std::endl;
//           ABORT("there");
//        }
//        ppos=epos;
//    }
    return Grad.size()?
                psiX.cwiseProduct(Grad.conjugate()).real()
              : res;
}

Eigen::MatrixXcd SchroedingerMachine::funcA(const Eigen::VectorXf &Elem, const Eigen::MatrixXcd &Psi,const Eigen::VectorXf Grad) const{
    Eigen::MatrixXd X(funcX(Elem));
    Eigen::MatrixXcd psi=Eigen::Map<const Eigen::MatrixXcd>(reinterpret_cast<const std::complex<double>*>(funcPhi(X,Psi,Eigen::MatrixXd()).data()),
                                                            X.rows(),X.cols());
    Eigen::MatrixXcd grad=Eigen::MatrixXcd::Random(X.rows(),X.cols());
    test(X,Psi,grad);
    return (psi.conjugate().cwiseProduct(psi));
}

Eigen::MatrixXcd SchroedingerMachine::funcB(const Eigen::VectorXf &Elem, const Eigen::MatrixXcd &Psi, const Eigen::VectorXf Grad) const {
    return funcA(Elem,Psi).cwiseProduct(funcW(Elem));
}

Eigen::MatrixXd SchroedingerMachine::LossL2(const Eigen::VectorXf &Elem, const Eigen::MatrixXcd &Psi, const Eigen::VectorXf Grad) const
{
    if(weigBasic().size()==0)setBasic(Elem.size());
    if(Grad.size()){
        Eigen::MatrixXd m(Elem.size(),1);
        //        m=funcB(Elem,Psi,Grad);
    }
    else{
        Eigen::MatrixXd m(1,1);
        m(0,0)=funcB(Elem,Psi).lpNorm<1>();
        return m;
    }
}

static Eigen::MatrixXd weigTrain;
static void setWeig(){
    if(weigTrain.size())return;
    weigTrain.resize(_trainBas.front()->order(),_trainBas.size());
    for(int m=0;m<_trainBas.size();m++){
        std::vector<double>node,weig;
        _trainBas[m]->dvrRule(node,weig);
        std::vector<double>v(_trainBas[m]->valNodes());
        if(m==0)v.insert(v.begin(),0.);
        if(m==_trainBas.size())v.push_back(0.);
        for(int i=0;i<weigTrain.rows();i++)weigTrain(i,m)=weig[i];
    }
}

static int cnt=0;
static Eigen::MatrixXcd grad;
static Eigen::MatrixXd dx;
static bool good=false;
void SchroedingerMachine::test(const Eigen::MatrixXd &X, const Eigen::MatrixXcd &Psi, const Eigen::MatrixXcd Grad) const{
    if(Grad.size()==0)return;
    if(cnt>100)ABORT(Sstr+"here"+cnt);
    cnt++;

//    if(not good)
    grad=Eigen::MatrixXcd::Random(Grad.rows(),Grad.cols());
    grad(0,0)=0.;
    for(int n=1;n<grad.cols();n++){
        grad(grad.rows()-1,n-1)=grad(0,n);
    }
    grad(grad.rows()-1,grad.cols()-1)=0.;

//    if(not good)
    dx=0.01*Eigen::MatrixXd::Random(X.rows(),X.cols());
    dx(0,0)=0.;
    for(int n=1;n<dx.cols();n++){
        dx(dx.rows()-1,n-1)=dx(0,n);
    }
    dx(dx.rows()-1,dx.cols()-1)=0;
    Eigen::MatrixXd d0(dx.cwiseProduct(funcPhi(X,Psi,grad)));
    double gana=d0.sum();

    Eigen::MatrixXd v0;
    v0=funcPhi(X+dx,Psi,Eigen::MatrixXcd())-funcPhi(X-dx,Psi,Eigen::MatrixXcd());
    v0=0.5*v0.cwiseProduct(Eigen::Map<const Eigen::MatrixXd>(reinterpret_cast<const double*>(grad.data()),grad.rows()*2,grad.cols()));
    double gnum=v0.sum();
    if(abs(gana/gnum-1.)>1.e-2){
        std::ofstream f0("f0");
        f0<<"# x,num,ana"<<std::endl;
         Eigen::MatrixXd vv=funcPhi(X,Psi,Eigen::MatrixXcd());
        Eigen::MatrixXd ww=funcPhi(X+dx,Psi,Eigen::MatrixXcd());
        std::cout<<"ana,num: "<<gana/gnum<<" "<<gana<<" "<<gnum<<std::endl;
        Eigen::MatrixXd res(dx.rows(),dx.cols());
        for(int j=0;j<dx.cols();j++)
            for(int i=0;i<dx.rows();i++){
                f0<<X(i,j)<<" "<<v0(2*i,j)+v0(2*i+1,j)<<" "<<d0(i,j)<<" "<<vv(2*i,j)+vv(2*i+1,j)<<" "<<ww(2*i,j)+ww(2*i+1,j)<<std::endl;

                Eigen::MatrixXd dy=Eigen::MatrixXd::Zero(dx.rows(),dx.cols());
                dy(i,j)=0.001;
                dy(0,0)=0.;
                dy(dy.rows()-1,dy.cols()-1)=0;
                gana=dy.cwiseProduct(funcPhi(X,Psi,grad)).sum();
                Eigen::MatrixXd v0=funcPhi(X,Psi,Eigen::MatrixXcd());
                dy+=X;
                Eigen::MatrixXd v1=funcPhi(dy,Psi,Eigen::MatrixXcd());
                gnum=(v1-v0).cwiseProduct(Eigen::Map<const Eigen::MatrixXd>(reinterpret_cast<const double*>(grad.data()),Grad.rows()*2,Grad.cols())).sum();
                res(i,j)=gnum/gana;
                if((gana==0.) != (gnum==0.))Sout+"not both zero"+gana+gnum+Sendl;
                if((gana!=0. and abs(gnum/gana-1.)<1.e-7) or (abs(gnum/gana-1.)>1.e-1 and abs(gnum)>1.e-6)){
                    Sout+"gnum==0"+i+j+gnum/gana+gnum+gana+Sendl;
//                    std::cout<<"v0"<<v0-v1<<std::endl;
                }
            }
//        std::cout<<res<<std::endl;
    }
    else good=true;
}
