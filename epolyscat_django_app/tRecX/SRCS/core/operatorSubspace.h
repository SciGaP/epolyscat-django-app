// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORSUBSPACE_H
#define OPERATORSUBSPACE_H

#include <memory>
#include "operatorTree.h"
#include "coefficients.h"

class ProjectSubspace;

/** \ingroup Structures */

///@brief Operator restricted to subspace
///
/// P Op P, where P is class ProjectSubspace
class OperatorSubspace: public OperatorTree{
    std::shared_ptr<const ProjectSubspace> _proj;
    const OperatorAbstract* _op;
    std::unique_ptr<Coefficients> _C,_D;
public:
    OperatorSubspace(const OperatorAbstract* Op); ///< infer projection from Op's Index (for Subspace&Complement hybrids)
    OperatorSubspace(const OperatorAbstract* Op, std::shared_ptr<const ProjectSubspace> Proj);
    std::shared_ptr<const ProjectSubspace> projector() const {return _proj;}
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
};


#endif // OPERATORSUBSPACE_H
