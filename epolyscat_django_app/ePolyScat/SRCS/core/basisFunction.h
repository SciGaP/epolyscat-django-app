// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef FUNCTIONSET_H
#define FUNCTIONSET_H

#include <vector>
#include "tools.h"
#include "constants.h"
#include "useMatrix.h"
#include "orthopol.h"

/** \ingroup Basissets */
/// Functions scaled to interval [0,1] or [0,infty) (used by many 1d BasisSet's)
class BasisFunction{
    friend class BasisSet;
    friend class BasisSetDef;
    friend class BasisFunctionRL;
    friend class BasisFunctionSqrt;
public:
    static void Test(bool Print); ///< perform basic tests
    void plot(std::string Dir) const; ///< create plot

public:
    virtual ~BasisFunction(){}
    BasisFunction():nameStr("undefined"),order(0),leftInfty(false),rightInfty(false){}
    BasisFunction(std::string NameStr,unsigned int Order);

    UseMatrix val(const UseMatrix & X) const; ///< return values on X grid
    UseMatrix der(const UseMatrix & X) const; ///< return derivatives on X grid

    /// base and weight for N-point quadrature on basis function
    virtual void quadRule(const unsigned int & N, UseMatrix & QuadX,UseMatrix & QuadW) const {ABORT("implement for "+nameStr);}

    /// DVR (=lobatto or radau) quadrature rule for BasisFunction
    virtual void dvrRule(UseMatrix & QuadX,UseMatrix & QuadW) const {ABORT("DVR rule not implemented)");}

    /// test various function kinds
    virtual void test(bool Print) const;

    /// true if function(infty)=0
    virtual bool asympZero() const{return false;}
    virtual bool leftZero() const{return false;}
    virtual bool rightZero() const{return false;}

    static bool asympZero(std::string Name);

    /// return string describing function kind
    virtual std::string name(bool Print=true) const {return nameStr;}
    unsigned int order; ///< number of functions

protected:
    const bool leftInfty;
    const bool rightInfty;
    bool checkStandardInterval(double X, bool ZeroOutside) const {
        //HACK
        return true;
        if(ZeroOutside) return true;
        if((not leftInfty and X<0.-1.e-14) or (not rightInfty and X>1.+1.e-14))
            ABORT("argument of function out of range, "+nameStr+": "+tools::str(X));
    }

    /// name string
    const std::string nameStr;

    /// parameters for function set (meaning is function dependent)
    std::vector<double> par;

    /// values and derivatives at a set of points X
    void valDer(const UseMatrix & X, UseMatrix & Val, UseMatrix & Der, bool ZeroOutside=false) const;

    /// values and derivatives at a set of points X
    virtual void valDer(const std::complex<double> & X, std::vector<std::complex<double> >  & Val, std::vector<std::complex<double> >  & Der, bool ZeroOutside=false) const
    {ABORT("implement valDer for "+nameStr);}

    /// table quad*[kind][nPoints][i] ... i'th quadrature point and weights in nPoints rule for kind
   static std::map<std::string,const BasisFunction*> tableNew;
    static const BasisFunction * get(std::string NameStr,unsigned int Order,std::vector<double>Par); ///< get basis function pointer from table (create if needed)
};

class BasisFunctionIndex:public BasisFunction{
    std::vector<std::vector<std::complex<double> > > vals;
public:
    BasisFunctionIndex(unsigned int Order);
    BasisFunctionIndex(const std::vector<int>&Vals);
protected:
    void valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const;
};

class BasisFunctionIncoming:public BasisFunction{
public:
    BasisFunctionIncoming(unsigned int Order):BasisFunction("incoming",Order){ABORT("re-implement \"incoming\"");}
};

class BasisFunctionGrid:public BasisFunction{
private:
    std::vector<double> points,weights;
public:
    BasisFunctionGrid(std::vector<double> Points, std::vector<double> Weights):BasisFunction("grid",Points.size()),points(Points),weights(Weights){}
protected:
    void valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const;
};

/// BasisFunction multiplied by R^l
class BasisFunctionRL: public BasisFunction{
private:
    const BasisFunction * basFun;
public:
    BasisFunctionRL(const BasisFunction* BasFun,std::vector<double>Par=std::vector<double>(0));
    bool rightZero() const {return basFun->rightZero();}
    bool leftZero()  const {return par[0]>0.;}
    bool asympZero() const {return basFun->asympZero();}
    void quadRule(const unsigned int &N, UseMatrix &QuadX, UseMatrix &QuadW) const {
        basFun->quadRule(N+(unsigned int)(par[0]),QuadX,QuadW);
    }
protected:
    void valDer(const UseMatrix &X, UseMatrix &Val, UseMatrix &Der, bool ZeroOutside) const;
};

/// BasisFunction multiplied by R^l
class BasisFunctionSqrt: public BasisFunction{
private:
    const BasisFunction * basFun;
public:
    BasisFunctionSqrt(const BasisFunction* BasFun);
    bool rightZero() const {return basFun->rightZero();}
    bool leftZero()  const {return false;}
    bool asympZero() const {return basFun->asympZero();}
    void quadRule(const unsigned int &N, UseMatrix &QuadX, UseMatrix &QuadW) const {
        basFun->quadRule(N+(unsigned int)(par[0]),QuadX,QuadW);
    }
    void dvrRule(UseMatrix &QuadX, UseMatrix &QuadW) const {basFun->dvrRule(QuadX,QuadW);}
protected:
    void valDer(const std::complex<double> & X, std::vector<std::complex<double> >  & Val, std::vector<std::complex<double> >  & Der, bool ZeroOutside=false) const;
    void valDer(const UseMatrix &X, UseMatrix &Val, UseMatrix &Der, bool ZeroOutside) const;
};


/// Cos and Sin
class BasisFunctionCosSin: public BasisFunction{
public:
    BasisFunctionCosSin(unsigned int Order,std::string Name="CosSin"):BasisFunction(Name,Order){
        for(unsigned int k=0;k<Order;k++){
            if(k%2==0)par.push_back(-double(    k/2));
            else      par.push_back( double((k+1)/2));
        }
    }
    bool rightZero() const {return false;}
    bool leftZero()  const {return false;}
    bool asympZero() const {return false;}
    void quadRule(const unsigned int &N, UseMatrix &QuadX, UseMatrix &QuadW) const {
        QuadX=UseMatrix(N,1);
        for(unsigned int k=0;k<N;k++)QuadX(k)=k/double(N);
        QuadW=UseMatrix::Constant(N,1,1/double(N));
    }
    void dvrRule(UseMatrix &QuadX, UseMatrix &QuadW) const {quadRule(order,QuadX,QuadW);}
protected:
    void valDer(const std::complex<double> &X, std::vector<std::complex<double> >&Val, std::vector<std::complex<double> >&Der, bool ZeroOutside) const;
};

/// complex, plane wave
class BasisFunctionPlaneWave: public BasisFunction {
    double deltaK;
public:
    BasisFunctionPlaneWave(unsigned int Order, std::vector<double> Par=std::vector<double>(),unsigned int Periods=1)
        :BasisFunction("PlaneWave",Order),deltaK(2*math::pi*Periods){
        if(Par.size()==0){
            for(unsigned int k=0;k<Order;k++){
                if(k%2==0)par.push_back(-double(    k/2));
                else      par.push_back( double((k+1)/2));
            }
        } else {
            if(Par.size()!=1)
                ABORT("specify first starting exponent, nothing else");
            par.clear();
            for(unsigned int k=0;k<Order;k++)
                par.push_back(Par[0]+k);
        }

    }
    bool rightZero() const {return false;}
    bool leftZero()  const {return false;}
    bool asympZero() const {return false;}
    void quadRule(const unsigned int &N, UseMatrix &QuadX, UseMatrix &QuadW) const {
        QuadX=UseMatrix(N,1);
        for(unsigned int k=0;k<N;k++)QuadX(k)=k/double(N);
        QuadW=UseMatrix::Constant(N,1,1/double(N));
    }
    void dvrRule(UseMatrix &QuadX, UseMatrix &QuadW) const {quadRule(order,QuadX,QuadW);}
protected:
    void valDer(const std::complex<double> &X, std::vector<std::complex<double> >&Val, std::vector<std::complex<double> >&Der, bool ZeroOutside) const;
};

/// OrthogonalPolynomial x sqrt(weight)
class BasisFunctionOpolWeig: public BasisFunction{
private:
    const OrthogonalPolynomial * oPol;
public:
    BasisFunctionOpolWeig(const OrthogonalPolynomial * OPol,unsigned int Order)
        :BasisFunction(OPol->name()+"*sqrt(w)",Order),oPol(OPol){}

    std::string name(bool Print=true) const {return oPol->name()+"*sqrt(w)";}
    bool rightZero() const {return oPol->upperBoundary()> DBL_MAX/2.;}
    bool leftZero()  const  {return oPol->lowerBoundary()<-DBL_MAX/2.;}
    bool asympZero() const {return rightZero() or leftZero();}
    void quadRule(const unsigned int &N, UseMatrix &QuadX, UseMatrix &QuadW) const;
    void dvrRule(UseMatrix &QuadX, UseMatrix &QuadW) const;
protected:
    void dvrRule(std::vector<double> &Node, std::vector<double>  &Weig) const;
    void valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const;
};

/// DVR type of basis function
class BasisFunctionDvr: public BasisFunctionOpolWeig{
    //-------------------------------------------------
    // WARNING:
    // DVR functions are NOT in ascending DVR point order
    // this is for running through the general FEM setup
    // where the first two functions are assumed to be
    // suitable to give lower/upper margin values
    //-------------------------------------------------
    std::vector<double> dvrNodes,qNorms;
public:
    BasisFunctionDvr(const OrthogonalPolynomial * OPol,unsigned int Order);
    std::string name(bool Print=true) const {return "DVR:"+BasisFunctionOpolWeig::name();}
protected:
    void valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const;
};

#endif // FUNCTIONSET_H
