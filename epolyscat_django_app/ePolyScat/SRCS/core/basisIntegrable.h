// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef BASISINTEGRABLE_H
#define BASISINTEGRABLE_H

#include <memory>

#include "basisAbstract.h"

#include "jacobian.h"
#include "complexScaling.h"

class UseMatrix;

/// an integrable basis must have an arbitrary point quadrature
class BasisIntegrable: public BasisAbstract
{
protected:
    double _lowBound,_upBound; ///< integration range
    std::shared_ptr<const Jacobian> _jac;
    std::shared_ptr<const ComplexScaling> _comSca;

public:
    static double infty;
    virtual ~BasisIntegrable(){}
    BasisIntegrable();
    BasisIntegrable(std::string Def);
    BasisIntegrable(double LowBound, double UpBound,std::string Jacob="1", ComplexScaling ComSca=ComplexScaling())
        :_lowBound(LowBound),_upBound(UpBound),
          _comSca(std::shared_ptr<const ComplexScaling>(new ComplexScaling(ComSca)))
    {
        if(Jacob=="1")_jac=std::shared_ptr<Jacobian1>(new Jacobian1(0.));
        else DEVABORT("awaits implementation");
    }

    virtual unsigned int lowerMargin() const {return 0;} ///< first function is on left boundary (should go to BasisIntegrable)
    virtual unsigned int upperMargin() const {return size()-1;}///< last function is on right boundary (should go to BasisIntegrable)

    virtual std::string str(int Level=0) const;
    virtual std::string strDefinition() const;
    virtual unsigned int order() const=0; ///< quadrature order for exact evaluation of overlap
    virtual void quadRule(int N, UseMatrix & QuadX, UseMatrix & QuadW) const;

    virtual void quadRule(int N, std::vector<double> &QuadX, std::vector<double> &QuadW) const=0;
    std::vector<std::complex<double> > val(double X) const;
    std::vector<std::complex<double> > der(double X) const;

    /// values and derivatives at points X, ZeroOutside: set = 0, where X outside range of basis function
    ///
    /// sorting: X-index runs faster, interpreted as column-major matrix: mat(i,n)= bas[n](X[i])
    virtual void valDer(const std::vector<std::complex<double> > & X,
                        std::vector<std::complex<double> > & Val,
                        std::vector<std::complex<double> > & Der, bool ZeroOutside=false) const=0;
    /// alternate call with strictly real coordinate poitns
    void valDerD(const std::vector<double > & X,
                 std::vector<std::complex<double> > & Val,
                 std::vector<std::complex<double> > & Der, bool ZeroOutside=false) const;


    UseMatrix val(const std::vector<double> &Grid) const; ///< val(i,k): value f_k(coordinate[i])
    UseMatrix val(const UseMatrix & Coordinates, bool ZeroOutside=false) const; ///< val(i,k): value f_k(coordinate[i])
    UseMatrix der(const UseMatrix & Coordinates, bool ZeroOutside=false) const; ///< derivatives at a coordinate grid

    /// old style version of valDer
    virtual void valDer(const UseMatrix & X,UseMatrix & Val, UseMatrix & Der, bool ZeroOutside=false) const;

    virtual double upBound() const {return _upBound;}       ///< return upper boundary of interval
    virtual double lowBound() const{return _lowBound;}      ///< return lower boundary of interval
    virtual std::vector<double> intervals() const {return {_lowBound,_upBound};} ///< points where basis is not analytic
    const std::shared_ptr<const Jacobian> jacobian() const {return _jac;}
    const std::shared_ptr<const ComplexScaling> complexScaling() const {return _comSca;}

    std::complex<double> eta() const;// {return _comSca->etaX(0.5*(upBound()+lowBound()));}
    bool isAbsorptive() const {return eta().imag()!=0.;}
    void plot(std::string Dir) const;
    std::vector<std::vector<std::complex<double> > > transZeroValDer(double Q) const; ///< unitary transformation to valDer(Q)=0
};

#endif // BASISINTEGRABLE_H
