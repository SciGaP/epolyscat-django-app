// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include <stdio.h>      /* printf */
#include <stdlib.h>
#include <complex>      // abs
#include <algorithm>    // std::max
#include <string>
#include <iostream>
#include <vector>

#include "str.h"

#include "orthopol.h"
#include "orthogonalDerived.h"

#include "useMatrix.h"
#include "tMath.h"
#include "constants.h"
#include "abort.h"
#include "tools.h"

// library returns quadarature rules
#include "qtAlglib.h"
#include "qtAlglib.h"

#include "integrate.h"
#include "vectorReal.h"


using namespace alglib;
using namespace std; 

shared_ptr<OrthogonalPolynomial> OrthogonalPolynomial::copyFactory(const OrthogonalPolynomial *Other){
    if(0!=dynamic_cast<const OrthogonalLaguerre*>(Other))
        return shared_ptr<OrthogonalPolynomial>(
                    new OrthogonalLaguerre(dynamic_cast<const OrthogonalLaguerre*>(Other)->alfa()));
    else if(0!=dynamic_cast<const OrthogonalLegendre*>(Other))
        return shared_ptr<OrthogonalPolynomial>(new OrthogonalLegendre());
    else if(0!=dynamic_cast<const OrthogonalChebychev*>(Other))
        return shared_ptr<OrthogonalPolynomial>(new OrthogonalChebychev());
    else
        ABORT("not in copyFactory: "+Other->name());
    return 0;
}


vector<double> OrthogonalPolynomial::val(int N, const double X) const{
    vector<double> V,dum;
    valDer(N,X,V,dum);
    return V;
}
vector<double> OrthogonalPolynomial::der(int N, const double X) const{
    vector<double> D,dum;
    valDer(N,X,dum,D);
    return D;
}

void OrthogonalPolynomial::valDer(int N, const double X, std::vector<double> &V, std::vector<double> &D) const {
    vector<long double> v,d;

    v.resize(N);
    d.resize(N);
    if(N>0){
        startValues(X,v[0],d[0]);
    }
    if (N>1) {
        v[1]=(a(1) + b(1)* X) * v[0];
        d[1]=(a(1) + b(1)* X) * d[0] + b(1) * v[0];
    }
    if (N>2)
    {
        for (int n=2; n<N;n++){
            v[n] = (a(n) +  b(n)*X)*v[n-1]               + c(n)*v[n-2];
            d[n] = (a(n) +  b(n)*X)*d[n-1] + b(n)*v[n-1] + c(n)*d[n-2];
        }
    };
    D.clear();
    V.clear();
    for(unsigned int n=0;n<v.size();n++){
        V.push_back(v[n]);
        D.push_back(d[n]);
    }
}
void OrthogonalPolynomial::alphaBetaMu(unsigned int N, vector<double> & Alpha, vector<double> &Beta, double & Mu) const{
    Alpha.clear();
    Beta.clear();
    if(N>0){
        Alpha.push_back(-a(1)/b(1));
        Beta.push_back(0.);
    }
    for(unsigned int n=2;n<N+1;n++){
        Alpha.push_back(-a(n)/b(n));
        Beta.push_back(-c(n)/(b(n)*b(n-1)));
    }
    Mu=normsq(0);
}

void OrthogonalPolynomial::quadrature(int N, vector<double>& X, vector<double>& W) const {

    // NOTE: alglib uses its own internal data structures

    alglib::ae_int_t n=N+addQuad();
    alglib::ae_int_t info;
    alglib::real_1d_array xq;
    alglib::real_1d_array wq;
    if(name().find("AssocLegendre")==0)
        alglib::gqgenerategausslegendre(n,info,xq,wq);
    else{
        alglib::real_1d_array alpha;
        alglib::real_1d_array beta;
        vector<double> stdAlpha,stdBeta;
        double mu0;
        alphaBetaMu(n,stdAlpha,stdBeta,mu0);
        alpha.setcontent(stdAlpha.size(),stdAlpha.data());
        beta.setcontent(stdBeta.size(),stdBeta.data());
        alglib::gqgeneraterec(alpha,beta,mu0,n,info,xq,wq);
        switch(info){
        case -3:ABORT("internal eigenproblem solver hasn't converged"+name());
        case -2:ABORT(Str("Beta[i]<=0:")+stdBeta);
        case -1:ABORT(Str("incorrect N was passed: ","")+n);
        }
    }

    // transfer to standard container
    X.resize(N+addQuad());
    W.resize(N+addQuad());
    for (size_t i=0; i<W.size(); i++) {
        X[i]=xq[i];
        if(lowerBoundary()!=-DBL_MAX and upperBoundary()!=DBL_MAX){
            if(abs(X[i]-lowerBoundary())<(upperBoundary()-lowerBoundary())*1.e-10)X[i]=lowerBoundary();
            if(abs(X[i]-upperBoundary())<(upperBoundary()-lowerBoundary())*1.e-10)X[i]=upperBoundary();
        }
        W[i]=wq[i];
    };
}
void OrthogonalPolynomial::quadratureKind(std::string Kind, int N, vector<double>& X, vector<double>& W) const {

    // NOTE: alglib uses its own internal data structures
    alglib::real_1d_array alpha;
    alglib::real_1d_array beta;
    vector<double> stdAlpha,stdBeta;
    double mu0;
    alglib::ae_int_t n=N;
    alglib::ae_int_t info=0;
    alglib::real_1d_array xq;
    alglib::real_1d_array wq;

    // for lobatto or radau
    alphaBetaMu(n,stdAlpha,stdBeta,mu0);
    alpha.setcontent(stdAlpha.size(),stdAlpha.data());
    beta.setcontent(stdBeta.size(),stdBeta.data());

    //std::cout << (Str("kind= ")+Kind +" b= " + tools::str(((OrthogonalCustom*)this)->_b)) << std::endl;
    
    if     (Kind=="Gauss")  alglib::gqgeneraterec(alpha,beta,mu0,                              n,info,xq,wq);
    else if(Kind=="RadauLeft")  alglib::gqgenerategaussradaurec(alpha,beta,mu0,lowerBoundary(),n,info,xq,wq);
    else if(Kind=="RadauRight") alglib::gqgenerategaussradaurec(alpha,beta,mu0,upperBoundary(),n,info,xq,wq);
    else if(Kind=="Lobatto")alglib::gqgenerategausslobattorec(alpha,beta,mu0,lowerBoundary(),upperBoundary(),n,info,xq,wq);
    else ABORT("undefined quadrature kind="+Kind);

    switch(info){
    case -3:ABORT(Str("internal eigenproblem solver hasn't converged")+name()+Kind+upperBoundary()+lowerBoundary());
    case -2:ABORT("Beta[i]<=0: "+tools::str(stdBeta,2,","));
    case -1:ABORT("quadrature rules need >=3 points, have: "+tools::str(N));
    }

    // transfer to standard container
    X.resize(xq.length());
    W.resize(wq.length());
    for (size_t i=0; i<W.size(); i++) {
        X[i]=xq[i];
        W[i]=wq[i];
    };
}

void OrthogonalPolynomial::quadratureWithEnds(int N, vector<double>& X, vector<double>& W) const {

    // NOTE: alglib uses its own internal data structures
    alglib::real_1d_array alpha;
    alglib::real_1d_array beta;
    vector<double> stdAlpha,stdBeta;
    double mu0;
    alglib::ae_int_t n=N;
    alglib::ae_int_t info = 0;
    alglib::real_1d_array xq;
    alglib::real_1d_array wq;

    // for lobatto or radau
    alphaBetaMu(n,stdAlpha,stdBeta,mu0);
    alpha.setcontent(stdAlpha.size(),stdAlpha.data());
    beta.setcontent(stdBeta.size(),stdBeta.data());

    if(lowerBoundary()<-DBL_MAX/2)
        ABORT("not implemented for lower boundary = Infty")
    else if(upperBoundary()>DBL_MAX/2)
        alglib::gqgenerategaussradaurec(alpha,beta,mu0,0.,n,info,xq,wq);
    else
        alglib::gqgenerategausslobattorec(alpha,beta,mu0,lowerBoundary(),upperBoundary(),n,info,xq,wq);

    switch(info){
    case -3:ABORT("internal eigenproblem solver hasn't converged (with ends)");
    case -2:ABORT("Beta[i]<=0: "+tools::str(stdBeta,2,","));
    case -1:ABORT("incorrect N was passed: "+tools::str(N));
    }

    // transfer to standard container
    X.resize(xq.length());
    W.resize(wq.length());
    for (size_t i=0; i<W.size(); i++) {
        X[i]=xq[i];
        W[i]=wq[i];
    };
}



// === end of functions defined for the abstract base class ===============
// NOTE: usually, code for the specific functions would follow, but we put this (inline) into the header
// exept for the quadrature rules

OrthogonalNassocLegendre::OrthogonalNassocLegendre(int M):OrthogonalPolynomial("AssocLegendre(normalized)"){
    m = std::abs(M);
    negative = M<0 and m%2!=0;
}
void OrthogonalNassocLegendre::startValues(long double X, long double &V, long double &D) const{
    // (see notes)
    V=std::sqrt(0.5*(long double)(2*m+1));
    if((m%2!=0) ^ negative)V=-V;
    for (unsigned int k=1;k<m+1;k++)V*=(2*k-1)/std::sqrt((2*k-1)*(2*k));
    V*=std::pow(std::sqrt(1.-X*X), m);

    // avoid floating dived by 0 at boundaries (exact value would be better)
    long double x=X;
    if(x<-1.+1.e-14)x=-1.+1.e-14;
    if(x> 1.-1.e-14)x=1.-1.e-14;
    D=-V*m*x/(1.-x*x);
}

long double OrthogonalNassocLegendre::b(int i) const {
    unsigned int l=i+m;
    return  std::sqrt((long double)((4*l*l-1))
                      /(long double)(l*l-m*m));
}
long double OrthogonalNassocLegendre::c(int i) const {
    unsigned int l=i+m;
    return -std::sqrt( (long double)((2*l+1)*((l-1)*(l-1)-m*m))
                       /(long double)((2*l-3)*(l*l-m*m)));
}

OrthogonalAssocLegendre::OrthogonalAssocLegendre(unsigned int M)
    :OrthogonalPolynomial("AssocLegendre"){negative=M<0;if(negative)m=-M;else m=M;}
double OrthogonalAssocLegendre::normsq(int I) const{
    double a=2./(long double)(2.*(I+m)+1);
    for(unsigned int k=I+1;k<I+2*m+1;k++)a*=k;
    return a;
}

void OrthogonalAssocLegendre::startValues(long double X, long double &V, long double &D) const{
    // P^m_m(x)=(-1)^m (2m-1)!! (1-x^2)^(m/2)
    if(negative){ABORT("negative m not implemented");}
    V=tMath::doubleFactorial( 2*m-1)*std::pow(std::sqrt(1.-X*X), m);
    if(m%2!=0)V=-V;
    D=-V*m*X/(1.-X*X);
}

bool OrthogonalPolynomial::test(int N, bool Print) {

    if(N>99){
        cout<<"OrthogonalPolynomial test limited to N<100\n";
        exit(1);
    };

    // get the N-point quadrature for these polynomials
    vector<double> x,w,v;
    quadrature(N,x,w);

    // compute the overlap matrix S
    vector<vector<double> > S(N,vector<double>(N,0.));
    for (int i=0;i<(int) x.size();i++) {
        v=val(N,x[i]);
        for (int m=0;m<(int)v.size(); m++) {
            for (int n=0;n< (int)v.size(); n++) {
                S[m][n]+=v[m]*v[n]*w[i];
            };
        };
    };

    // check orthogonality and print norms
    vector<long double> qnorm;
    for (unsigned int n=0;n<v.size();n++)qnorm.push_back(1./std::sqrt(normsq(n)));
    long double err=0.,nerr=0.;
    for (int m=0; m< (int) v.size(); m++) {
        for (int n=0; n< (int) v.size(); n++) {
            if(m==n) {
                nerr=max(nerr,(abs(S[m][n]-normsq(n))*(qnorm[n]*qnorm[m])));
            } else {
                err=max(err,abs(S[m][n])*qnorm[n]*qnorm[m]);
            };
        };
    };

    if (max(nerr,err)>tolerance()) {
        if(Print) cout<<"SERIOUS ERROR: "+name()+" N="<<N<<", Maximal deviation from orthogonality, norm: "<<setw(10)<<err<<" "<<nerr<<endl;
        return false;
    } else {
        if(Print)cout<<"OK "+name()+" N="<<N<<", max error: "<<setw(10)<<err<<endl;
        return true;
    };
}

void OrthogonalPolynomial::Test(bool Print) {
    cout<<"\nTest orthogonality and norm"<<endl;

    OrthogonalLegendre().test(50,Print);
    OrthogonalJacobi(-0.4,-0.4).test(50,true);
    OrthogonalJacobi(0.,-0.4).test(50,true);
    OrthogonalJacobi(0.,0.1).test(50,true);
    OrthogonalJacobi(0.,2.).test(50,true);
    OrthogonalLaguerre().test(50,Print);
    OrthogonalLaguerre(1).test(3,Print);
    OrthogonalLaguerre(3).test(50,Print);
    OrthogonalLaguerre(21).test(50,Print);
    OrthogonalChebychev().test(10,true);

//    OrthogonalDerived(5,1.,5.,"x",lin).test(5,true);
//    OrthogonalDerived(5,1.,5.,"legendre",one).test(5,true);
//    OrthogonalDerived(5,1.,5.,"assLeg(1)",asLeg1).test(5,true);
//    OrthogonalDerived(5,1.,5.,"assLeg(4)",asLeg4).test(5,true);

    OrthogonalNassocLegendre(0).test(5,Print);
    OrthogonalNassocLegendre(-1).test(5,Print);
    OrthogonalNassocLegendre(-5).test(5,Print);
}


long double OrthogonalJacobi::a(int i) const {if(A*A-B*B==0) return 0.;   return (2*i+A+B-1)*(A*A-B*B)/(2*i*(i+A+B)*(2*i+A+B-2));}
long double OrthogonalJacobi::b(int i) const {if(A+B==0)return double(2*i-1)/double(i); return (2*i+A+B-1)*(2*i+A+B)/(2*i*(i+A+B));}
long double OrthogonalJacobi::c(int i) const {if(A==0 and B==0)return -double(2*i-1)/double(i); return -2*(i+A-1)*(i+B-1)*(2*i+A+B)/(2*i*(i+A+B)*(2*i+A+B-2));}
double OrthogonalJacobi::weight(double X) const {return pow(1-X,A)*std::pow(1+X,B);}
double OrthogonalJacobi::derWeight(double X) const {
    if(A==0.)return  B*std::pow(1+X,B-1.);
    if(B==0.)return -A*std::pow(1-X,A-1.);
    return -A*pow(1-X,A-1.)*std::pow(1+X,B)+B*pow(1-X,A)*std::pow(1+X,B-1.);
}


/**
 * orthogonal polynomials for the measure x-b on [-1;+1]
 * note: b must not be the zero of any legendre-polynomial
 * (which is automatically the case if b less than -1)
 * 
 * recurrence relation: Legendre
 * (k+1)*Pk+1 = (2k+1)*X*Pk - k * Pk-1
 * (k+2)*Pk+2 = (2k+3)*X*Pk+1 - (k+1)*Pk
 *
 * define Bk = Pk(b)
 * 
 * the custom polynomials are
 * Qk := 1/(X-b)*(Pk+1 - Bk+1/Bk * Pk)
 * 
 * Qk+1 - (2k+3)/(k+2) * X*Qk
 * =1/(X-b)/(k+2)*[(k+2)*Pk+2-(k+2)*Bk+2/Bk+1*Pk+1 
 *  - (2k+3)*X*Pk+1 + (2k+3)*X*Bk+1/Bk*Pk ]
 * =1/(X-b)/(k+2)*[-(k+1)*Pk - (k+2)*Bk+2/Bk+1*Pk+1
 *  + (2k+3)*X*Bk+1/Bk* Pk]
 * = a* [-(k+2)*Bk+2/Bk+1 * Pk+1 + (2k+3)*Bk+1/Bk * X*Pk - (k+1)*Pk]
 * 
 * a = 1/(k+2) * 1/(X-b)
 * c = -(k+2)*Bk+2/Bk+1
 * d = (2k+3)*Bk+1/Bk
 * 
 * ... = a*c*[Pk+1] + a*d * [X*Pk] - a*(k+1)*Pk
 *     = a*c*[Pk+1 - Bk+1/Bk * Pk]
 *      +a*c*[Bk+1/Bk * Pk]
 *      +a*d/(2k+1)*[(2k+1)*X*Pk - k*Pk-1]
 *      +a*d/(2k+1)*[k*Pk-1]
 *      -a*(k+1)*Pk
 *     = a*c*Qk + a*c*[Bk+1/Bk * Pk]
 *      +a*d/(2k+1)*(k+1)*Pk+1
 *      +a*d/(2k+1)*k* Pk-1
 *      -a*(k+1)*Pk
 * 
 * e = a*d/(2k+1)*(k+1)
 * 
 * ... = a*c*Qk + e*[Pk+1-Bk+1/Bk*Pk] + e*[Bk+1/Bk*Pk]
 *      +a*c*[Bk+1/Bk * Pk] - a*(k+1)*Pk
 *      +a*d/(2k+1)*k* Pk-1
 *     = (a*c+e)*Qk
 *     + [e*Bk+1/Bk + a*c*Bk+1/Bk - a*(k+1)]*Pk
 *     + a*d/(2k+1)*k* Pk-1
 *     = (a*c+e)*Qk
 *     + [e*Bk+1/Bk + a*c*Bk+1/Bk - a*(k+1)]*[Pk-Bk/Bk-1*Pk-1]
 *     + { ... } // should be 0
 * 
 * Qk+1 = {(2k+3)/(k+2)} * X*Qk
 *      + {-Bk+2/Bk+1
 *         +(k+1)/(k+2)*(2k+3)/(2k+1) * Bk+1/Bk } Qk
 *      + { (k+1)/(k+2)*(2k+3)/(2k+1) * Bk+1*Bk+1/Bk/Bk
 *         -Bk+2/Bk+1
 *         -(k+1)/(k+2) } Qk-1
 * 
 * These factors are b(k+1), a(k+1), c(k+1)
 * 
 * start values: 
 * Q-1 = 0
 * Q0 = 1
 * Q1 = X + 1/(3*b) = X*Q0 + 1/(3*b) * Q0
 * 
 * ?? find fast expression for Tk or the terms in Tk (perhaps doesn't exist)
 */
OrthogonalCustom::OrthogonalCustom(double __a, double __b): OrthogonalPolynomial("custom"), _a(__a), _b(__b) {
    /*OrthogonalLegendre P = OrthogonalLegendre();
    vector<double> B = P.val(10+2, _b);
    for(int k=2; k<10; k++) {
        std::cout << k << ": " << B[k] << " " << a(k) << " " << b(k) << " " << c(k) << std::endl;
    }*/
}

double OrthogonalCustom::normsq(int I) const { // \int P_I(x) P_I(x) w(x) dx
    /*
     |Qn|² = \int (Pn+1 - Bn+1/Bn * Pn)² 1/(x-b) dx 
     *     = -Bn+1/Bn * \int Pn * (Pn+1-Bn+1/Bn*Pn) 1/(x-b) dx  // because Pn+1 orthogonal to polynomials of lower degree
     *     = -Bn+1/Bn * \int Pn * ((2k+1)/(k+1)*X*Pn + ...) * 1/(x-b) dx
     *     = -Bn+1/Bn * \int Pn*Pn * (2k+1)/(k+1) dx +
     *       -Bn+1/Bn * \int Pn*(...)/(x-b) dx
     * (...)/(x-b) is polynomial of degree less than n
     *     = -Bn+1/Bn*(2k+1)/(k+1) * \int Pn*Pn dx
     */
    OrthogonalLegendre P = OrthogonalLegendre();
    vector<double> B = P.val(I+2, _b);
    return _a*(-B[I+1]/B[I] * (2.*I+1.)/(I+1.) * P.normsq(I));
}
long double OrthogonalCustom::lowerBoundary() const {return -1;}
long double OrthogonalCustom::upperBoundary() const {return +1;}
double OrthogonalCustom::weight(double X) const {
    return _a*(X-_b);
}
double OrthogonalCustom::derWeight(double X) const {
    return _a*1.0;
}

/*void OrthogonalCustom::startValues(long double X, long double & V, long double & D) const {
    V=1.,D=0;  ///< value and derivative of lowest degree function
}*/
    

long double OrthogonalCustom::a(int k) const {
    k -= 1;
    OrthogonalLegendre P = OrthogonalLegendre();
    vector<double> B = P.val(k+3, _b);
    if(k >= 0) {
//        return 1 - (k+1.)/(k+2.)*(2.*k+3.)/(2.*k+1.);
        return -B[k+2]/B[k+1] + (k+1.)/(k+2.)*(2.*k+3.)/(2.*k+1.) * B[k+1]/B[k];
    }
    else {
        DEVABORT("a(0) should not be needed in OrthogonalPolynomial::a");
        return 0;
    }
}
long double OrthogonalCustom::b(int k) const {
    k -= 1;
    return (2.*k+3.)/(k+2.);
}
long double OrthogonalCustom::c(int k) const {
    k -= 1;
    OrthogonalLegendre P = OrthogonalLegendre();
    vector<double> B = P.val(k+3, _b);
    if(k >= 1){
        return 1./(k+2.)*(2.*k+3.)*B[k+1]*B[k+1]/B[k]/B[k] / (2.*k+1)*(k+1) - B[k+2]/B[k] - (k+1.)/(k+2.);
    }
    else {
        DEVABORT("c(0),c(1) should not be needed in OrthogonalPolynomial::c");
        return 0;
    }
}

/*
template<class T>
OrthogonalScaled<T>::OrthogonalScaled(double a, T P): _a(a), _P(P) {}
template<class T>
double OrthogonalScaled<T>::normsq(int I) const {return _a * _P.normsq(I);}
template<class T>
long double OrthogonalScaled<T>::lowerBoundary() const {return _P.lowerBoundary();}
template<class T>
long double OrthogonalScaled<T>::upperBoundary() const {return _P.upperBoundary();}
template<class T>
double OrthogonalScaled<T>::weight(double X) const {return _a*_P.weight(X);}
template<class T>
double OrthogonalScaled<T>::derWeight(double X) const {return _a*_P.derWeight(X);}
template<class T>
long double OrthogonalScaled<T>::a(int i) const {return _P.a(i);}
template<class T>
long double OrthogonalScaled<T>::b(int i) const {return _P.b(i);}
template<class T>
long double OrthogonalScaled<T>::c(int i) const {return _P.c(i);}
*/
