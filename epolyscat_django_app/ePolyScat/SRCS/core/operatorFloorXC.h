// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORFLOORXC_H
#define OPERATORFLOORXC_H

#include "operatorFloor.h"
#include "operatorZG.h"
#include "operatorZD.h"


class OperatorTree;
class BasisIntegrable;
class BasisOrbitalNumerical;

/** \ingroup OperatorFloors */

///@brief exchange mean field operators for a given PotentialMultipole
///
/// mapping two radial intervals s in [a,b] -> r in [c,d] the general form
///
///      X[i,j](r,s) chi(s) =  V(r,s) rho[i,j](r,s)
///
/// where V(r,s) are potentials in a suitable quadrature representation and
///
///     rho[i,j](r,s)=int ds1..dsN Psi[i](r,s1,...,sN)Psi[j](s,s1,...,sN)
///
/// is the sub-block of the generalized 1-particle reduced density matrix
/// the density matrix is represented as
///
///      rho[i,j]=sum[lmkn] Ylm(r) conj(Ykn(s)) radial[lmkn](r,s)
///
class OperatorFloorXC : public OperatorZG
{
    Eigen::MatrixXcd _mat;
public:
    static void postProcess(OperatorTree * Op, std::string Pot, const std::vector<std::vector<Eigen::MatrixXcd> > &RhoIJ,
                            const BasisOrbitalNumerical *IOrb);

    OperatorFloorXC(std::string Pot, const Index *AIndex, const Index *BIndex, std::complex<double> Multiplier);
    bool hasBeenSetUp() const {return _mat.size()>0;}
    Eigen::MatrixXcd & mat(){return _mat;}
    static void test(const OperatorTree* Op, const BasisOrbitalNumerical &IOrb);
private:
    void finalize(const OperatorTree *OpLeaf);

};

#endif // OPERATORFLOORXC_H
