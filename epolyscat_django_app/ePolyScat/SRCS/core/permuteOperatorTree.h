// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PERMUTEOPERATORTREE_H
#define PERMUTEOPERATORTREE_H

#include <complex>
#include <string>

#include "operatorTree.h"
#include "tree.h"
#include "index.h"
#include "coefficients.h"
#include "operatorFloor.h"

//#define _PERMUTE_OPERATOR_TREE_NEWCODE_

/** \ingroup Structures */
/// @brief Create OperatorTree associated to permuted index structures
///
/// Given an OperatorTree and a permutation of the iIndex and jIndex structures, this class
/// creates a permuted operator tree which has as i/jIndex the permuted indices. Functionality
/// associated to permuting Coefficients is also encapsulated here.
/// 
/// Warning! Has never reached a well-tested "doubt-free" stage. Specifically when usingthe NEWCODE
/// flag, be sure to check, whether permutation worked. Handling of timeDepFac might be an issue in
/// certain cases. See comments in cpp
class PermuteOperatorTree
{
    friend class TimePropagatorOutput;
    // Unpermuted objects
    const Index* iIndex;
    const Index* jIndex;
    const OperatorTree* optree;

    // Permuted objects
    const Index* iPermutedIndex;
    const Index* jPermutedIndex;
    OperatorTree* permutedOptree;
    unsigned int permutedFloorLevel;

    // Permutations of iIndex & jIndex
    std::vector<unsigned int> iPermutation;
    std::vector<unsigned int> jPermutation;


    /*
     * Helpers
     */
private:
    class OperatorFloorSingle: public OperatorFloor{
        std::complex<double> val;

    protected:
        void axpy(std::complex<double> alpha, const std::vector<std::complex<double> > &x,
                  std::complex<double> beta, std::vector<std::complex<double> > &y) const{
            y[0] = alpha*val*x[0]+beta*y[0];
        }
        void axpy(const std::complex<double> & alpha, const std::complex<double>* x, unsigned int SizX,
                  const std::complex<double> & beta,        std::complex<double>* y, unsigned int SizY) const{
            y[0] = alpha*val*x[0]+beta*y[0];
        }
    public:
        OperatorFloorSingle(): OperatorFloor(1,1,"ZS"){}
        OperatorFloorSingle(std::complex<double> _val): OperatorFloor(1,1,"ZS"), val(_val) {}


        void pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const{}
    };

    static void move(std::vector<unsigned int>& perm, unsigned int from, unsigned int to);
    static std::vector<unsigned int> invertPermutation(std::vector<unsigned int> permutation);

    static void permuteCoeffs(const std::vector<unsigned int> permutation, const Index* srcIndex, const Coefficients& src, const Index* targetIndex, Coefficients& target);


	static std::vector<unsigned int> classifyStructureForLeaf(const OperatorTree* optree);
	static std::vector<unsigned int> classifyStructure(const OperatorTree* optree);
#ifdef _PERMUTE_OPERATOR_TREE_NEWCODE_
    static OperatorTree* buildOptreeFromMatrix(
			const UseMatrix& matrix, std::string name, const Index* iIndex, const Index* jIndex, const Index* iIndexRoot, const Index* jIndexRoot);
    static OperatorTree* alignToIndex(const OperatorTree* optree, const Index* iIndex, const Index* jIndex);
	static void insertFloorLevel(OperatorTree* optree);

#else
    static OperatorTree* buildOptreeFromMatrix(
        const UseMatrix& matrix, std::string name, const Index* iIndex, const Index* jIndex, const Index* iIndexRoot, const Index* jIndexRoot, std::complex<double>* timeDepFac);
    static void findIndexChildren(OperatorTree* optree, const Index* iIndex, const Index* jIndex, std::vector<OperatorTree*>& result);
    static void findFloorChildren(OperatorTree* optree, const Index* iIndex, const Index* jIndex, std::vector<OperatorTree*>& result);
	static OperatorTree* alignToIndex(const OperatorTree *optree, const Index *iIndex, const Index *jIndex);
    void placeInPermutedOptree0(OperatorTree* permutedOptree0, const OperatorTree* leaf, OperatorFloor* floor);

#endif

public:
    PermuteOperatorTree(const OperatorTree* _optree); ///< Standard constructor

    void iMove(unsigned int from, unsigned int to); ///< In iIndex move the from'th level to the to'th level
    void jMove(unsigned int from, unsigned int to); ///< In jIndex move the from'th level to the to'th level
    void move(unsigned int from, unsigned int to);  ///< Move in both iIndex and jIndex

    void permute(const Index* iPermutedIndex=0, const Index* jPermutedIndex=0); ///< Create permuted operator tree

    OperatorTree* getPermutedOperatorTree() const{ return permutedOptree; } ///< Getter

    /// Permute lhs coefficients
    void iPermuteCoeffs(const Coefficients& src, Coefficients& target) const{
        permuteCoeffs(iPermutation, iIndex, src, iPermutedIndex, target);
    }

    /// Permute rhs coefficients
    void jPermuteCoeffs(const Coefficients& src, Coefficients& target) const{
        permuteCoeffs(jPermutation, jIndex, src, jPermutedIndex, target);
    }

    /// Unpermute lhs coefficients, that is, given permuted coefficients, returns original ones.
    void iUnpermuteCoeffs(const Coefficients& src, Coefficients& target) const{
        permuteCoeffs(invertPermutation(iPermutation), iPermutedIndex, src, iIndex, target);
    }

    /// Unpermute rhs coefficients, that is, given permuted coefficients, returns original ones.
    void jUnpermuteCoeffs(const Coefficients& src, Coefficients& target) const{
        permuteCoeffs(invertPermutation(jPermutation), jPermutedIndex, src, jIndex, target);
    }

    /// Set the permuted floor level
    void setFloorLevel(unsigned int depth){
        permutedFloorLevel = depth;
    }

    /// Check if multplication with permuted operator yields the same result as multiplication with unpermuted one.
    /// Will ABORT if not the case.
    void check();

    /// @brief Simpple test of the class
    static void test();
};


#endif // OPERATORTREE_H
