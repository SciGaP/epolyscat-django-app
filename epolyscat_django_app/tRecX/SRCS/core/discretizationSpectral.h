// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONSPECTRAL_H
#define DISCRETIZATIONSPECTRAL_H

#include <memory>
#include <string>
#include "discretizationDerived.h"
#include "operatorTree.h"
#include "constrainedView.h"

class EigenSolver;
class OperatorDiagonal;
class ProjectSubspace;

/** \ingroup Discretizations */

/// eigenvectors of operator
class DiscretizationSpectral: public DiscretizationDerived{
    friend class DiscretizationSpectralProduct;
    std::string constString;

    int selectNmin(const Index* Idx) const{return 0;}
    int selectNmax(const Index* Idx) const;

//    void indexAndMaps(EigenSolver *Slv, Index * &Idx, OperatorTree *&FromParent, OperatorTree *&ToParent);
    OperatorTree * _mapConstructor(bool From, const Discretization *Parent, const Index* SIndex, const Index *EIndex,
                                   const std::vector<Coefficients*> & Evec);

    std::unique_ptr<ProjectSubspace> _project;

protected:
    OperatorDiagonal * spectralValues;
    std::string _selectionCriterion;
public:
    virtual ~DiscretizationSpectral();
    const std::vector<std::complex<double>>& eigenvalues() const;// {return spectralValues->diagonal();}
    OperatorDiagonal * spectralOper() const {return spectralValues;}
    DiscretizationSpectral(const Discretization *D, std::string Criterion);

    /// \brief construct spectral discretization from constrained view on an operator (OBSOLESCENT)
    DiscretizationSpectral(const Discretization *D,
                           const OperatorAbstract* Op,
                           double Emin=-DBL_MAX,
                           double Emax=DBL_MAX,
                           bool excludeEnergyRange = false
            );
    /// \brief spectral discretization for constrained view on an operator
    DiscretizationSpectral(const OperatorAbstract* Op,
                           double Emin=-DBL_MAX,
                           double Emax=DBL_MAX,
                           bool excludeEnergyRange = false
            );

    void check(const OperatorAbstract* Op) const;
    const OperatorAbstract* projector() const {return _project.get();}

};

#endif // DISCRETIZATIONSPECTRAL_H
