// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "quadratureRule.h"

#include "qtEigenDense.h"
#include "integrate.h"
#include "tools.h"
#include "basisIntegrable.h"
#include "algebra.h"
#include "basisSetDef.h"
#include "basisMat1D.h"

using namespace std;

namespace QuadratureRule {

static const BasisIntegrable* integrandBasis=0;
static const Algebra* integrandAlgebra=0;


static Eigen::MatrixXcd integrandBasisAlgebra(const vector<double>&X)
{
    if(X.size()!=1)ABORT("need vector length =1");
    vector<complex<double> > val,dum;
    integrandBasis->valDer(vector<complex<double> >(1,X[0]),val,dum);

    // conjg(bas[i](q)) * algebra(Q)
    return Eigen::Map<Eigen::MatrixXcd>(val.data(),1,val.size()).adjoint()*integrandAlgebra->val(X[0]);
}

static Eigen::MatrixXcd integrandBasisBasis(const vector<double>&X)
{
    if(X.size()!=1)ABORT("need vector length =1");
    vector<complex<double> > val,dum;

    integrandBasis->valDer(vector<complex<double> >(1,X[0]),val,dum);


    // outer product of values
    Eigen::MatrixXcd res(val.size(),val.size()*2);
    res.block(0,0,val.size(),val.size())
            =Eigen::Map<Eigen::MatrixXcd>(val.data(),1,val.size()).adjoint()*Eigen::Map<Eigen::MatrixXcd>(val.data(),1,val.size());

    // X * (product of values)
    res.block(0,val.size(),val.size(),val.size())=X[0]*res.block(0,0,val.size(),val.size());

    return res;
}

class Int:public Integrate::Tools
{
public:
    Int( double AccRel=1.e-12, double AccAbs=1.e-12,
         std::vector<std::vector<unsigned int> > NQuad=std::vector<std::vector<unsigned int> >(0),
         std::string Kind="GaussLegendre",
         std::string KindInf="GaussLaguerre")
        : Integrate::Tools(AccRel,AccAbs,NQuad,Kind,KindInf){}

    Eigen::MatrixXcd nDim(const std::vector<std::vector<double > > Vol,
                          const std::function<Eigen::MatrixXcd(const std::vector<double> &)> Func,
                          const std::vector<double> Params=std::vector<double>(0)
            ) {return Integrate::NDim<Eigen::MatrixXcd,double,Int>(*this,Vol,Func,Params);}

    Eigen::MatrixXcd recursive(const std::vector<std::vector<double> > Vol,
                               const std::function<Eigen::MatrixXcd(const std::vector<double> &)> Func,
                               const std::vector<double> Params=std::vector<double>(0)
            ) {return Integrate::Recursive<Eigen::MatrixXcd,double,Int>(*this,Vol,Func,Params);}
};

vector<complex<double> > integralsBasisAlgebra(const BasisIntegrable * Bas, const Algebra * Alg){
    integrandBasis=Bas;
    integrandAlgebra=Alg;

    // compute exact overlap <Bas|Algebra>
    Int integ; // create the integration object

    // 1-d integral over first and possibly further sub-intervals
    Eigen::MatrixXcd res=integ.recursive({{Bas->intervals()[0],Bas->intervals()[1],1.}},integrandBasisAlgebra);
    for(int k=2;k<Bas->intervals().size();k++)res+=integ.recursive({{Bas->intervals()[k-1],Bas->intervals()[k],1.}},integrandBasisAlgebra);

    return vector<complex<double> >(res.data(),res.data()+res.size());
}


vector<complex<double> > integralsBasisBasis(const BasisIntegrable * Bas){
    integrandBasis=Bas;

    // compute exact overlap <Bas|Algebra>
    Int integ; // create the integration object

    // 1-d integral over first and possibly further sub-intervals
    Eigen::MatrixXcd res=integ.recursive({{Bas->intervals()[0],Bas->intervals()[1],1.}},integrandBasisBasis);
    for(int k=2;k<Bas->intervals().size();k++)res+=integ.recursive({{Bas->intervals()[k-1],Bas->intervals()[k],1.}},integrandBasisBasis);

    return vector<complex<double> >(res.data(),res.data()+res.size());
}


void pointsAndWeights(const BasisIntegrable * Bas, std::vector<double> & Points, std::vector<double> & Weights)
{
    // WARNING: does not work on general (in particular translationally invariant orthogonal functions)

    integrandBasis=Bas;

    // compute exact overlap <Bas|Bas> and multiplication <Bas|Q|Bas>
    Int integ; // create the integration object

    // 1-d integral over first and possibly further sub-intervals
    Eigen::MatrixXcd res=integ.recursive({{Bas->intervals()[0],Bas->intervals()[1]}},integrandBasisBasis);
    for(int k=2;k<Bas->intervals().size();k++)res+=integ.recursive({{Bas->intervals()[k-1],Bas->intervals()[k]}},integrandBasisBasis);

    // diagonalize - eigenvalues of Q are quadrature points
    int siz=Bas->size();

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> slv(res.block(0,siz,siz,siz),res.block(0,0,siz,siz));

    // quadrature point q_k is the k'th eigenvalue, weight w_k=1/|eigenfunction[k](q_k)|^2
    Points.clear();
    Weights.clear();
    for(int k=0;k<siz;k++){
        vector<complex<double> > valQ,dum;
        Points.push_back(slv.eigenvalues()(k));

        Bas->valDer({Points.back()},valQ,dum);
        complex<double> valQk=0.;
        for(int l=0;l<siz;l++)valQk+=valQ[l]*slv.eigenvectors()(l,k);
        Weights.push_back(1./std::norm(valQk));
    }
}

static int testcnt=0;
void test(){
    if(testcnt++>0)return;
    std::vector<BasisSetDef> defs;
    defs.push_back(BasisSetDef(10,3.,2.,"polynomial",false,true,false,Coordinate::fromString("Rn")));
    defs.push_back(BasisSetDef(10,0.,10.,"besselCoulomb",false,false,false,Coordinate::fromString("Rn"),
    {0,-1},ComplexScaling(),false,{0}));

    for(BasisSetDef basDef: defs){
        const BasisIntegrable * bas=dynamic_cast<const BasisIntegrable*>(BasisAbstract::factory(basDef));
        if(bas==0)ABORT("fail");
        vector<double>pts,wgs;
        pointsAndWeights(bas,pts,wgs);
        vector<complex<double> > ovrX=integralsBasisBasis(bas);
        vector<complex<double> > ovrQ(ovrX.size(),0.);
        for(int k=0;k<pts.size();k++){
            vector<complex<double> >val,dum;
            bas->valDer({pts[k]},val,dum);
            for(int j=0,ij=0;j<bas->size();j++){
                complex<double> ValWeig=val[j]*wgs[k];
                for(int i=0;i<bas->size();i++,ij++){
                    ovrQ[ij]+=conj(val[i])*ValWeig;
                }
            }
        }
        UseMatrix ox=UseMatrix::UseMap(ovrX.data(),bas->size(),bas->size());
        UseMatrix oq=UseMatrix::UseMap(ovrQ.data(),bas->size(),bas->size());
        oq-=ox;
        if(oq.maxAbsVal()>1.e-12*ox.maxAbsVal())Sstr+oq.str("diff",0)+Sendl;
    }
    ABORT("stop tests");

}

}
