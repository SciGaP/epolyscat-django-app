// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORHARTREE_H
#define OPERATORHARTREE_H


#include "operatorFloor.h"
#include "operatorZD.h"

class OperatorTree;
class BasisOrbitalNumerical;

/** \ingroup OperatorFloors */

///@brief Hartree potential for channels {i,j}
///
///      D{i,j} chi(r) =  int ds V(r,s) rho{i,j}(r,s) chi[r]
///
/// where V(r,s) are potentials in a suitable quadrature representation (class MultipolePotential) and
///
///     rho{i,j}(r,s)=int ds1..dsN Psi[i](r,s1,...,sN)Psi[j](s,s1,...,sN)
///
/// is the sub-block of the generalized 1-particle reduced density matrix
/// the density matrix is represented as (class DensityMatrix1)
///
///      rho{i,j}(r,s)=sum[la,ma,lb,mb] Y[la,ma](r) conj(Y[lb,mb](s)) radial[la,ma,lb,mb](r,s)
///
/// Op.def():
/// "<Hartree[PotName]>", "<Hartree>" defaults to PotName="CoulombEE"
///
/// "<HartreeRelative{c}>": for this the c'th channel diagonal potential is subtracted from all channel diagonal pots
///
/// "<Hartree@OrbDef@DensMat>": OrbDef...valid orbital definition, DensMat...density matrix name, input by 'Matrix: ..."
///
class OperatorHartree : public OperatorZD
{
    friend class OperatorFloorHF;
    Eigen::MatrixXcd _mat; // temporary storage for accumulating floor matrix (resized by LeafC)
    void finalize(const OperatorTree *OpLeaf);
    static void moMatrixElements(const OperatorTree *Op, const BasisOrbitalNumerical *IOrb, const Eigen::MatrixXcd & RhoIJ); // compare results of application against direct matrix elements
public:

    /// calculate the hartree potential floors; call after complete OperaterHartree has been set up (efficiency)
    static void postProcess(OperatorTree * Op, const std::vector<std::vector<Eigen::MatrixXcd> > &RhoIJ,
                            const BasisOrbitalNumerical *IOrb);

    /// single hartree potential floor (needs postProcess)
    OperatorHartree(std::string PotChan, const Index *IIndex, const Index *JIndex,std::complex<double> Multiplier);

    Eigen::MatrixXcd & mat(){return _mat;} ///< access to temporary matrix storage
    bool hasBeenSetUp() const {return _mat.size()>0 or dat!=0;} ///< skip in PatchRadial (because already set up)

    static void test(const OperatorTree* Op, const BasisOrbitalNumerical &IOrb);
};

#endif // OPERATORHARTREE_H
