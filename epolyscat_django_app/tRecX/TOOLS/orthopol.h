// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ORTHOPOL_H
#define ORTHOPOL_H

#include "abort.h"
#include <stdio.h>      /* printf */
#include <stdlib.h>
#include <complex>      // abs
#include <algorithm>    // std::max
#include <string>
#include <iostream>
#include <vector>
#include <iomanip>
#include <memory>
#include <float.h>

// library returns quadarature rules
#include "qtAlglib.h"
#include "qtAlglib.h"

#include "integrate.h"

//WARNING: tgmath seems to lead to a range of conflicts
//#include "tgmath.h"

/** \addtogroup Functions
 * @{
 */
/// \brief base class for orthogonal polynomial
class OrthogonalPolynomial{
    static unsigned int currentOrder;
    static const OrthogonalPolynomial * currentPol;

    /// alpha and beta for standard recurrence (see turff.pdf)
    void alphaBetaMu(unsigned int N, std::vector<double> &Alpha, std::vector<double> &Beta, double & Mu) const;

    const std::string thisName;
protected:

    /// short notation for orthogonal function kinds
    virtual unsigned int addQuad() const {return 0;}
    virtual void startValues(long double X, long double & V, long double & D) const {V=1.,D=0;} ///< value and derivative of lowest degree function
    virtual double tolerance() const {return 1.e-13;} ///< tolerance for orthonormality tests, some do not pass at this level

    // switch between different quadrature kinds
    void quadratureKind(std::string Kind, int N, std::vector<double>& X, std::vector<double>& W) const;
    /// functions for recurrence relations
    virtual long double a(int i) const =0;
    virtual long double b(int i) const =0;
    virtual long double c(int i) const =0;

    OrthogonalPolynomial():thisName("undefined"){}
    OrthogonalPolynomial(std::string Name):thisName(Name){}

public:
    std::string name() const {return thisName;}

    virtual double weight(double X) const=0;    ///< weight function
    virtual double derWeight(double X) const=0; ///< derivative of weight function

    virtual double normsq(int I) const =0;
    virtual long double lowerBoundary() const {return -DBL_MAX;}
    virtual long double upperBoundary() const {return  DBL_MAX;}

    /// return values/derivatives up to degree N-1 (order N)
    virtual void valDer(int N, const double X, std::vector<double>& V, std::vector<double>& D) const;

    std::vector<double> val(int N, const double X) const;
    std::vector<double> der(int N, const double X) const;

    void quadrature(int N, std::vector<double>& X, std::vector<double>& W) const;
    void quadratureWithEnds(int N, std::vector<double>& X, std::vector<double>& W) const; ///< N-point quadratur rule containing the end points (unless +- infty)

    void quadratureGauss(int N, std::vector<double>& X, std::vector<double>& W) const{quadratureKind("Gauss",N,X,W);}
    void quadratureRadauLeft( int N, std::vector<double>& X, std::vector<double>& W) const{quadratureKind("RadauLeft",N,X,W);}
    void quadratureRadauRight(int N, std::vector<double>& X, std::vector<double>& W) const{quadratureKind("RadauRight",N,X,W);}
    void quadratureLobatto(int N, std::vector<double>& X, std::vector<double>& W) const{quadratureKind("Lobatto",N,X,W);}

    /// use quadrature rules to test up to degree N-1 (order N)
    bool test(int N, bool Print=false);
    static void Test(bool Print);

    /// copy OrthogonalPolynomial of suitable type
    static std::shared_ptr<OrthogonalPolynomial> copyFactory(const OrthogonalPolynomial* Other);
};

// ==== end of abstract base class header, below there are specifications ====

class OrthogonalLegendre:public OrthogonalPolynomial{
public:
    OrthogonalLegendre():OrthogonalPolynomial("legendre"){}
    inline double normsq(int I) const {return 2./(2*I+1);}
    inline long double lowerBoundary() const {return -1.;}
    inline long double upperBoundary() const {return  1.;}
    inline double weight(double X) const {return 1.;}
    inline double derWeight(double X) const {return 0.;}
private:
    // defined here in place, actually will be put directly into the code at compile time
    inline long double a(int i) const {return 0.;}
    inline long double b(int i) const {return (2.*double(i)-1.)/double(i);}
    inline long double c(int i) const {return   (-double(i)+1.)/double(i);}
};

/**
 * monic version Qn of OrthogonalLegendre Pn
 * Qn = (n!)^2 * 2^n / (2n)! * Pn
 */
class OrthogonalLegendreMonic:public OrthogonalPolynomial{
public:
    OrthogonalLegendreMonic():OrthogonalPolynomial("legendreMonic"){}
    inline double normsq(int I) const {
        return 2./(2*I+1) * pow(tgamma(I+1),4)*pow(2,2*I)/pow(tgamma(2*I+1),2);
    }
    inline long double lowerBoundary() const {return -1.;}
    inline long double upperBoundary() const {return  1.;}
    inline double weight(double X) const {return 1.;}
    inline double derWeight(double X) const {return 0.;}
private:
    // defined here in place, actually will be put directly into the code at compile time
    inline long double a(int i) const {return 0.;}
    inline long double b(int i) const {return 1.;}
    inline long double c(int i) const {return -(i-1.)*(i-1.)/(4.*(i-1.)*(i-1.)-1.);}
};

/**
 * measure w(x)=a*(x-b) on [-1;1]
 * 
 * Pk(x) = Monic Legendre-Polynomials
 * 
 * Qk(x) = 1/(x-b) * [Pk+1(x) - Pk+1(b)/Pk(b) * Pk(x)] = 
 *   orthogonal polynomials for the measure (x-b) and thus also for a*(x-b)
 */
class OrthogonalCustom:public OrthogonalPolynomial{
public:
    double _a,_b;
public:
    OrthogonalCustom(double a, double b);
    double normsq(int I) const;
    long double lowerBoundary() const;
    long double upperBoundary() const;
    double weight(double X) const;
    double derWeight(double X) const;
private:
    long double a(int i) const;
    long double b(int i) const;
    long double c(int i) const;
};

class OrthogonalChebychev:public OrthogonalPolynomial{
public:
    OrthogonalChebychev():OrthogonalPolynomial("chebychev"){}
    inline double normsq(int I) const {if(I==0)return M_PI;return M_PI/2;}
    inline long double lowerBoundary() const {return -1.;}
    inline long double upperBoundary() const {return  1.;}
    inline double weight(double X) const {return 1./sqrt(abs(1.-X*X));}
    inline double derWeight(double X) const {return X/(pow(sqrt(abs(1.-X*X)),3));}
private:
    // defined here in place, actually will be put directly into the code at compile time
    inline long double a(int i) const {return 0;}
    inline long double b(int i) const {return 2;}
    inline long double c(int i) const {return -1;}
};

/// normalized associated Legendre polynomials
class OrthogonalNassocLegendre:public OrthogonalPolynomial{
public:
    OrthogonalNassocLegendre(int M);
    double normsq(int I) const {return 1.;}
    inline long double lowerBoundary() const {return -1.;}
    inline long double upperBoundary() const {return  1.;}
    inline double weight(double X) const {return 1.;} // note: historically, we have the defined the weight into the functions
    inline double derWeight(double X) const {return 0.;}
private:
    unsigned int m;
    bool negative;

    // defined here in place, actually will be put directly into the code at compile time
    inline long double a(int i) const {return 0.;}
    long double b(int i) const;
    long double c(int i) const;

    unsigned int addQuad() const {return 0;}
    void startValues(long double X, long double &V, long double &D) const;
};

class OrthogonalLaguerre:public OrthogonalPolynomial{
private:
    int _alfa;
    inline long double a(int i) const {return (2.*double(i)-1.+_alfa)/double(i);}
    inline long double b(int i) const {return	                 -1./double(i);}
    inline long double c(int i) const {return   (-double(i)+1.-_alfa)/double(i);}
public:
    OrthogonalLaguerre(int Alfa=0):OrthogonalPolynomial("laguerre("+tools::str(Alfa)+")"),_alfa(Alfa){}
    inline double normsq(int I) const {
        double res=1.;
        for(int n=I+1;n<I+_alfa+1;n++)res*=double(n);
        return res;
    }
    inline long double lowerBoundary() const {return 0.;}
    inline double weight(double X) const {return exp(-X)*pow(X,_alfa);}
    inline double derWeight(double X) const {
        if(_alfa==0)return -exp(-X);
        if(X==0.)return 0.;
        return exp(-X)*pow(X,_alfa-1)*(_alfa-X);
    }
    double alfa() const {return _alfa;}
};

class OrthogonalAssocLegendre:public OrthogonalPolynomial{
public:
    OrthogonalAssocLegendre(unsigned int M);
    double normsq(int I) const;
    inline long double lowerBoundary() const {return -1.;}
    inline long double upperBoundary() const {return  1.;}
    inline double weight(double X) const {return 1.;} // note: historically, we have the defined the weight into the functions
private:
    unsigned int m;
    bool negative;

    // defined here in place, actually will be put directly into the code at compile time
    inline long double a(int i) const {return 0.;}
    inline long double b(int i) const {return (2.*(long double)(i+m)-1.)/(long double)(i);}
    inline long double c(int i) const {return (-(long double)(i+2*m)+1.)/(long double)(i);}

    unsigned int addQuad() const {return 0;}
    void startValues(long double X, long double &V, long double &D) const;
};

class OrthogonalHermite:public OrthogonalPolynomial{
public:
    OrthogonalHermite();
private:
    inline long double a(int i) const {return  0.;          }
    inline long double b(int i) const {return  2.;          }
    inline long double c(int i) const {return -2.*double(i-1);}
    inline double weight(double X) const {return exp(-X*X);}
};

class OrthogonalJacobi:public OrthogonalPolynomial{
    const double A,B;
    const double gA,gB,gAB;
    double gammaAB(int I) const {
        double gn=(gA*gB)/gAB;
        for(int n=1;n<=I;n++)gn*=(n+A)*(n+B)/(n*(n+A+B));
        return gn;
    }
public:
    OrthogonalJacobi(double Alfa, double Beta)
        :OrthogonalPolynomial("jacobi["+tools::str(Alfa,2)+","+tools::str(Beta,2)+"]"),
          A(Alfa),B(Beta),gA(tgamma(A+1)),gB(tgamma(B+1)),gAB(tgamma(A+B+1)){}
    inline double normsq(int I) const {return std::pow(2.,A+B+1)/(2*I+A+B+1)*gammaAB(I);}
    inline long double lowerBoundary() const {return -1.;}
    inline long double upperBoundary() const {return 1.;}
    double weight(double X) const;// {return pow(1-X,A)*std::pow(1+X,B);}
    double derWeight(double X) const;// {return -A*pow(1-X,A-1.)*std::pow(1+X,B)+B*pow(1-X,A)*std::pow(1+X,B-1.);}
private:
    long double a(int i) const;
    long double b(int i) const;
    long double c(int i) const;
};
/** @} */
///
#endif
