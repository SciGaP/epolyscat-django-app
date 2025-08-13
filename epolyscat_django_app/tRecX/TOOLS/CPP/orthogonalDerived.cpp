// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "orthogonalDerived.h"

#include "algebra.h"
#include "useMatrix.h"

using namespace std;

// for testing generic
//static double asLeg4(double x){return pow(abs(1.-x*x),3);}
//static double asLeg1(double x){return abs(1.-x*x);}
//static double one(double x){return 1.;}
//static double lin(double x){return x;}
//static double sqx(double x){return x*x;}
//static double che(double x){return 1./std::sqrt(1-x*x);}

double OrthogonalDerived::weight(double X) const {return _aWeig->val(X).real()*_base->weight(X);}
double OrthogonalDerived::derWeight(double X) const {return _aDer->val(X).real()*_base->derWeight(X);}

OrthogonalDerived::OrthogonalDerived(unsigned int MaxDegree, const OrthogonalPolynomial *Base, string Weig, string Der)
    :OrthogonalPolynomial(Base->name()+"["+Weig+"]")
{
    _base=OrthogonalPolynomial::copyFactory(Base);
    _aWeig=shared_ptr<Algebra>(new Algebra(Weig));
    if(Der!="")_aDer=shared_ptr<Algebra>(new Algebra(Der));

    if(not _aWeig->isAlgebra())ABORT("weight function string is not Algebra: "+Weig);
    if(Der!="" and not _aDer->isAlgebra())ABORT("derivative function string is not Algebra: "+Der);

    vector<vector<double> > ovr=overlap(MaxDegree+1);
    construct(ovr);

    // compute norms (we might, just as well, renormalize to 1)
    _normsq.assign(MaxDegree,0.);
    std::vector<double> x,w;
    quadratureGauss(MaxDegree,x,w);
    for (size_t k=0;k<x.size();k++){
        std::vector<double> val,dum;
        valDer(MaxDegree,x[k],val,dum);
        for(size_t n=0;n<val.size();n++)_normsq[n]+=val[n]*val[n]*w[k];
    }
}

void OrthogonalDerived::construct(std::vector<std::vector<double> > Ovr){
    // where is the docu for this algorithm? (but seems to be working!)

    // get transformation to orthonormal polynomials
    vector<vector<double> > trans;
    tools::gramSchmidtTrans(Ovr,trans);

    // evaluate on a quadrature grid for the basis polynomials
    vector<double> q,w;
    _base->quadrature(Ovr.size(),q,w);

    // evaluate new polynomials at quadrature points
    vector<std::complex<double> > pq(Ovr.size()*q.size());
    for(unsigned int k=0;k<q.size();k++){
        vector<double> bq=_base->val(Ovr.size(),q[k]);

        // apply transformation and add into final values
        for(unsigned int n=0;n<Ovr.size();n++){
            for(unsigned int l=0;l<trans[n].size();l++)
                pq[n*Ovr.size()+k]+=trans[n][l]*bq[l];
        }
    }

    // calculate recurrence coefficients as D = Q X Q^1
    UseMatrix qix;
    qix=UseMatrix::UseMap(pq.data(),trans.size(),trans.size());
    UseMatrix ddd(qix);
    for(unsigned int k=0;k<q.size();k++)ddd.row(k)*=q[k];
    qix.solve(ddd);

    // transform to standard form of coefficients
    recCoe.assign(3,vector<double>(Ovr.size()-1));
    for(unsigned int n=0;n<Ovr.size()-1;n++){
        recCoe[0][n]       =  ddd(n,n).real()/ddd(n+1,n).real();
        recCoe[1][n]       =  -1.            /ddd(n+1,n).real();
        if(n>0)recCoe[2][n]=ddd(n-1,n).real()/ddd(n+1,n).real();

        // remove near-zeros
        double cMax=max(max(abs(recCoe[0][n]),abs(recCoe[1][n])),abs(recCoe[2][n]));
        for(unsigned int l=0;l<3;l++){
            if(abs(recCoe[l][n])<cMax*1.e-12)recCoe[l][n]=0.;
        }
        // make leading coefficient > 0
        if(recCoe[2][n]>0. or recCoe[0][n]<=0.)
            for(unsigned int l=0;l<3;l++)recCoe[l][n]=-recCoe[l][n];
    }


}

void OrthogonalDerived::verify(){
    OrthogonalLaguerre lag2(3);
    lag2.test(20,true);

    OrthogonalLaguerre lag0;
    OrthogonalDerived der2(10,&lag0,"Q*Q*Q","3*Q*Q");
    der2.test(10,true);

    vector<double> val0,der0,val1,der1;
    bool isOK;
    std::vector<double> rat;
    for(double x=der2.lowerBoundary();x<=min(der2.upperBoundary(),der2.lowerBoundary()+2.);x+=0.3){
        isOK=true;
        lag2.valDer(6,x,val0,der0);
        der2.valDer(val0.size(),x,val1,der1);
        if(rat.size()==0){
            for(size_t k=0;k<val0.size();k++)rat.push_back(val0[k]/val1[k]);
        }
        for(size_t k=0;k<val0.size();k++){
            if(abs(val0[k]-val1[k]*rat[k])>1.e-10*sqrt(lag2.normsq(k)))isOK=false;
            if(abs(der0[k]-der1[k]*rat[k])>1.e-10*sqrt(lag2.normsq(k)))isOK=false;
        }
        if(not isOK){
            Sstr+"ERROR at x="+x+Sendl;
            Sstr+"explicit"+val0+Sendl;
            Sstr+" generic"+val1+Sendl;
            isOK=true;
        }
    }
    std::cout<<"OK OrthogonalDerived verified for Laguerre"<<std::endl;
}


vector<vector<double> > OrthogonalDerived::overlap(unsigned int Degree){
    vector<double> x,w;
    _base->quadrature(Degree+20,x,w);

    vector<vector<double> > res(Degree+1,vector<double>(Degree+1,0.));
    for(unsigned int k=0;k<x.size();k++){
        double fk=weight(x[k])/_base->weight(x[k]);
        vector<double> v=_base->val(Degree+1,x[k]);
        for(unsigned int i=0;i<res.size();i++){
            for(unsigned int j=0;j<res[i].size();j++){
                res[i][j]+=fk*w[k]*v[i]*v[j];
            }
        }
    }
    return res;
}
