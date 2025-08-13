// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPPERATORRALL_H
#define OPPERATORRALL_H

#include <memory>

#include "operatorTree.h"

class OperatorTree;

/// \brief rALL - recursive Approximation by Local Low Rank decomposition
class OperatorRALL : public OperatorTree
{
    struct _term{
        /// triplet of (u, si, v)
        std::unique_ptr<OperatorTree> u,v;
        std::unique_ptr<OperatorAbstract>si;
        std::unique_ptr<Coefficients> tmpR,tmpL;
        _term(OperatorAbstract *Si, OperatorTree* U, OperatorTree* V):u(U),v(V),si(Si){
            if(si){
                tmpR.reset(new Coefficients(si->jIndex));
                tmpL.reset(new Coefficients(si->iIndex));
            }
            else {
                if(v)tmpR.reset(new Coefficients(v->iIndex));
                if(u)tmpL.reset(new Coefficients(u->jIndex));
            }
        }
    };
    std::vector<_term> _terms;
    const size_t maxTerms;
public:

    /// construct from complex symmetric OperatorTree with blocks at Height above oFloor
    OperatorRALL(OperatorTree& O /** complex symmetric operator, destroyed during construction */,
                 std::vector<size_t> Height={0} /** use blocks at Height above floor, length determines depth of recursion  */,
                 size_t BlockBand=1, /** use block-band as first term, band-width 2*Band-1, none for Band=0 */
                 size_t MaxTerms=1 /** maximal number of expansion terms */);

    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const override;
};

#endif // OPPERATORRALL_H
