// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PROJECTSUBSPACE_H
#define PROJECTSUBSPACE_H

#include <memory>
#include "operatorAbstract.h"
#include "index.h"
#include "eigenSolver.h"

class OperatorTree;
/** \ingroup Structures */
///@brief Project onto a subspace
///
/// sequence of maps from full to expansion in terms of orthonormal subspace basis and back
class ProjectSubspace : public OperatorAbstract
{
    std::vector<int> _sorting;
    Index* _subspaceIndex; // potential memory leak - found no solution for handling this any better
    std::shared_ptr<OperatorTree> _mapFrom,_mapTo; // also used in DiscretizationSpectral
    std::unique_ptr<Coefficients> _subspaceC;
    std::shared_ptr<Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>> _lu;
    void _construct(std::vector<const Coefficients*> Vectors, std::vector<const Coefficients*> Duals, size_t BeginOrthonormal=0);
public:
    ProjectSubspace(EigenSolver & Slv);
    ProjectSubspace(std::vector<Coefficients*> Vectors, std::vector<Coefficients*> Duals,size_t BeginOrthonormal=0);
    ProjectSubspace(std::vector<const Coefficients*> Vectors, std::vector<const Coefficients*> Duals, size_t BeginOrthonormal=0);
    std::shared_ptr<OperatorTree> mapFrom() const { return _mapFrom;} ///< map from the reference discretization to vector
    std::shared_ptr<OperatorTree> mapTo() const { return _mapTo;}     ///< map to the reference discretization from vector
    const Index* subspaceIndex() const {return _subspaceIndex;}
    std::shared_ptr<const Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>,Eigen::COLAMDOrdering<int>>> lu() const {return _lu;}
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    /// apply S P S^-1: properly project onto dual space
    virtual void applyDual(std::complex<double> A, const Coefficients &Dual, std::complex<double> B, Coefficients &Y) const;
    int dim() const {return _subspaceIndex->size();}

    bool verify() const;
    const std::vector<int> & sorting() const {return _sorting;}

    /// return the orbitals corresponding to projector
    std::vector<Coefficients> orbitals(std::vector<int> Select={} /** list of coefficients, defaults to all */);
    virtual ~ProjectSubspace();
};

/// map to dual space rather than rhs space
class ProjectSubspaceToDual: public OperatorAbstract{
    std::shared_ptr<ProjectSubspace> _proj;
    mutable Coefficients _tmp;
public:
    ProjectSubspaceToDual(std::shared_ptr<ProjectSubspace> Proj)
        :OperatorAbstract(Proj->name,Proj->idx(),Proj->jdx()){_proj=Proj;_tmp.reset(idx());}
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
        _proj->apply(A,Vec,0,_tmp);
        _tmp.idx()->overlap()->apply(1,_tmp,B,Y);
    }

};

#endif // PROJECTSUBSPACE_H
