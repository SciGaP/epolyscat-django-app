// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATION_COUPLED_LM_H
#define DISCRETIZATION_COUPLED_LM_H

#include "discretizationDerived.h"
#include "basisAbstract.h"

#include "operatorTree.h"

#include <vector>

/**
 * Discretization of \f$L^2(\mathbb R^3)\otimes L^2(\mathbb R^3)\f$ in terms of coupled angular/magnetic momentum
 * L, M. This is modeled as a derived discretization from one using individual angular/magnetic momenta.
 *
 * In other words, mapFromParent is the matrix of Clebsch-Gordan coefficients
 * \f[
 *      <l_1l_2LM|l_1m_1l_2m_2>\delta_{\dots},
 * \f]
 *
 * where \f$\delta_{\dots}\f$ means diagonal in all other (i. e. radial) indices.
 */
class DiscretizationCoupledLM: public DiscretizationDerived{
public:
    class IndexedBasisAbstract: public BasisAbstract{
        std::vector<int> indices;

    public:
        IndexedBasisAbstract(std::string Name, std::vector<int> Indices):
            BasisAbstract(Name),
            indices(Indices){}

        BasisAbstract* remove(const std::vector<int>& RemoveK) const override{
            std::vector<int> newIndices(indices);

            std::vector<int> removeK(RemoveK);
            std::sort(removeK.begin(), removeK.end());

            for(int i=removeK.size() - 1; i>=0; i--){
                newIndices.erase(newIndices.begin() + removeK[i]);
            } 

            return new IndexedBasisAbstract(name(), newIndices);
        }
        unsigned int size() const override{ return indices.size(); }
        double physical(int Index) const override{ return indices[Index]; }
    };


private:
    class MapFromParent: public OperatorTree{
    public:
        MapFromParent(const Index* IIndex, const Index* JIndex, double* floor_factor);
    };

    class MapToParent: public OperatorTree{
    public: 
        MapToParent(const Index* IIndex, const Index* JIndex, double* floor_factor);;
    };

    Index* setupIndexStructure(int M1, int L1, int M2, int L2, std::vector<int> idx, const Index* model);
public:
    DiscretizationCoupledLM(const Discretization* Parent);
};






#endif // DISCRETIZATION_COUPLED_LM_H
