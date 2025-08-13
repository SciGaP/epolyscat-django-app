// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef BASPROD_H
#define BASPROD_H

#include "basisNdim.h"
#include "coordinateTrans.h"

/** \ingroup Basissets */
/// tree product basis
///
/// quadrature grid either native to product basis or for reference BasisNdim
///\n
class BasisProd : public BasisNdim
{
    std::unique_ptr<const Index> _idx;
    mutable std::unique_ptr<Index> _gdx;
    const BasisNdim * ref;
    static std::map<const Index*,std::map<const Index*,BasisProd*>> all;

    // tensor product Prod <- Prod (x) Factor[0] (x) ... (x) Factor.back()
    std::vector<std::complex<double> > tensorMultiply(const std::vector<std::complex<double> > & LFac, const std::vector<std::complex<double> > &RFac) const;
    // evaluate all basis functions in tree-product Idx at one coordinate point
    void valDerAll(const Index * Idx, std::vector<double> Point, std::vector<std::vector<std::complex<double> > > &ValDer) const;

    // it becomes obvious that these function pointers should be in BasisNdim and matching actual functions are not needed
    // this is all strictly coordinate dependent and not much dependent on basis
    CoordinateMap _fromNdim,_toNdim;/// Ndim may be Phi.Eta.R, Phi.Rho (to be extended)
    JacobianMap _jacToNdim; ///< Ndim may be Phi.Eta.R, Phi.Rho (to be extended)
    IntegrationFactor _intFactor;
    NablaFactor _nabFactor;

    mutable std::vector<std::vector<double>> _quadMutable;
    mutable std::vector<double> _weigMutable;

protected:
    std::vector<double> toCartesian(const std::vector<double>&CoorNdim) const {return _fromNdim(CoorNdim);}
    std::vector<double> fromCartesian(const std::vector<double>&CoorNdim) const {return _toNdim(CoorNdim);}
    double absFactor(const std::vector<double> & CoorNdim) const{return _intFactor(CoorNdim);}
    double nablaFactor(const std::vector<double>&CoorNdim,int I) const{return _nabFactor(CoorNdim,I);}
    Eigen::MatrixXd jacobianToNdim(const std::vector<double>&CoorNdim) const;

public:
    /// BasisNdim with quadrature rule from BasisNdim and values adjusted for quadrature rule, i.e. quadCoor()!=ndimCoor()
    BasisProd(const Index *RootIdx, const BasisNdim *Ref);
    /// BasisNdim with product quadrature rule for bases in tree and product values, i.e. quadCoor()==ndimCoor()
    BasisProd(const Index *RootIdx, unsigned int MinQuad=10);
    static const BasisProd* factory(const Index* Root, const Index * Ref); ///< use quadrature grid of reference
    std::vector<std::complex<double> > operator()(std::vector<double> CoorNdim) const;
};

#endif // BASPROD_H
