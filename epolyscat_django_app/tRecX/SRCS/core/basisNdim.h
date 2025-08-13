// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef BASISNDIM_H
#define BASISNDIM_H

#include "basisAbstract.h"
#include <map>
#include "useMatrix.h"
#include "complexScaling.h"
#include "vectorValuedFunction.h"

class CoordinateTrans;
class OperatorDefinition;
class Index;
class OperatorTree;

/** \ingroup Basissets */

/// Basis function depending on several arguments
///
/// function values and partial derivatives at a list of n-dimensional quadrature points and weights
/// is stored in the class
///
/// quadrature points and partial derivatives are wrt to the reference coordinates quadCoor()
/// e.g. spherical with origin (0,0,0) or Cartesian
///
///
class BasisNdim : public BasisAbstract, public VectorValuedFunction
{
    static std::map<const std::string,std::map<const Index*,std::map<const Index*,UseMatrix> > > opers;
    static std::map<std::string,const BasisAbstract*> bases;
    static std::map<std::string,const OperatorTree*> basicOp;

protected:
    std::string _ndimCoor; // coordinates wrt to basis
    std::string _quadCoor; // coordinates for _quadGrid
    std::vector<std::vector<double> > _quadGrid; // quadrature points in terms of the main coordinate system
    std::vector<double> _quadWeig; // quadrature weight wrt to main coordinate systems (and integration measure)
    std::vector<std::vector<std::vector<std::complex<double>>>> _valDer; //  value and derivative of all basis functions at all quadrature points

    std::vector<ComplexScaling> _comSca;

    // essential properties of any Ndim basis
public:
    using BasisAbstract::name;
    /// definition of coordinates: transform CoorNdim to Cartesian
    virtual std::vector<double> toCartesian(const std::vector<double>&CoorNdim) const=0;
    /// (back-)transform from Cartesian to CoorNdim
    virtual std::vector<double> fromCartesian(const std::vector<double>&Cartesian) const=0;
protected:
    /// Jacobian d(cartesian)/d(ndim coords) at point expressed in CoordNdim
    virtual Eigen::MatrixXd jacobianToNdim(const std::vector<double>&CoorNdim) const=0;
    /// extra integration factor at point, typically sqrt(Jacobi-det) (see tsurff.pdf for details)
    virtual double absFactor(const std::vector<double>&CoorNdim) const=0;
    /// partial derivatives of absFactor wrt ndim coords at point CoorNdim
    virtual double nablaFactor(const std::vector<double>&CoorNdim,int I) const=0;


    /// quadrature grid and weights for tree-product basis
    void productGridWeig(const Index * Idx, std::vector<std::vector<double> > & Grid, std::vector<double> & Weig, int MinQuad,
                         std::vector<double> GridK=std::vector<double>(), double WeigK=1.);
    void mainQuadValDer(); ///< move quadrature rule, values, and gradient to quadrature coordinates
    std::string mainCoor(ReadInput &Inp);
public:
    static std::vector<const BasisNdim*> allNdim(std::string Name); /// get all in group form table
    static const BasisAbstract* factory(std::string Name); /// get from table
    static void read(ReadInput & Inp); ///< put into table from input
    /// evaluate matrix at operator floor involving at least one basis NDim
    static void matrix(const std::string &Op, const Index* IIndex, const Index* JIndex, UseMatrix & Mat, std::complex<double> Multiplier);

    std::string ndimCoor() const {return _ndimCoor;}
    std::string quadCoor() const {return _quadCoor;} ///< the coordinates used for the quadrature (grid and weight)
    unsigned int dim() const {return std::count(_quadCoor.begin(),_quadCoor.end(),'.')+1;} ///< spatial dimension
    unsigned int size() const {return (_valDer.size()?_valDer[0].size():0);} ///< number of basis functions
    std::string str(int Level=0) const;

    /// quadrature points wrt main coordinate system
    virtual const std::vector<std::vector<double> > & quadGrid() const {return _quadGrid;}
    virtual size_t quadSize() const {return _quadGrid.size();}
    virtual const std::vector<double> & quadGrid(size_t Pt) const {return _quadGrid[Pt];}
    /// quadrature weight wrt to main coordinate systems (and integration measure)
    virtual const std::vector<double> & quadWeig() const {return _quadWeig;}
    virtual double quadWeig(size_t Pt) const {return _quadWeig[Pt];}
    /// valDers()[pt][ibas][k]: at pt'th point, and ibas'th basis funtion  value [k=0] and partial derivatives [k>0] wrt. quadCoor() coordintes
    virtual const std::vector<std::vector<std::vector<std::complex<double>>>> &valDers() const {return _valDer;}
    virtual std::complex<double> value(size_t Pt, size_t Ibas){return valDers()[Pt][Ibas][0];}
    virtual std::complex<double> valPartial(size_t Pt, size_t Ibas, size_t Kpartial /** 0=value, 1,2,..partial derivatives*/) const {return valDers()[Pt][Ibas][Kpartial];}

    const Index* rootIndex(const Index* Idx) const; ///< return index root matching BasisNdim, 0 if no match

    bool isOrthonormal() const {return false;}

    void test();

    /// for VectorValuedFunction
    std::vector<std::complex<double> > operator()(std::vector<double> CoorNdim) const=0;
    std::string coordinates() const{return ndimCoor();} ///< coordinate for evaluating at arbitrary point
    unsigned int length() const{return size();} ///< length of vector

    /// temporary, until moved to Index or elswhere
    static std::vector<std::complex<double>> IndexAt(const Index* Idx, const std::vector<double> Coors);

};


#endif // BASISNDIM_
