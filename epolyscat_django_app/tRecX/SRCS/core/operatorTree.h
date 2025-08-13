// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORTREE_H
#define OPERATORTREE_H

#include <complex>
#include <string>

#include "qtEigenSparse.h"
#include "tree.h"
#include "operatorAbstract.h"
#include "coefficientsLocal.h"
#include "operatorDefinition.h"

class OperatorDefinition;
class OperatorFloor;
class UseMatrix;
class TransformOperatorAbstract;
class OperatorTucker;
class HMatrix;
class SparseMatrix;

/** \ingroup Structures */
/// new main operator type base on template class Tree
class OperatorTree: public Tree<OperatorTree>, public OperatorAbstract
{
    friend class OperatorAbstract;
    friend class HMatrix;
    friend class Parallel;
    friend class PermuteOperatorTree;
    friend class OperatorTucker;
    friend class TransformOperatorAbstract; // for setting floors
    friend class OperatorSVD; // for setting floors
    friend class Index;
    friend class DerivativeBlock; //HACK for changing parallelization
    friend class ParallelOperator; // that is OK
    friend class BasisOrbital; //HACK - should be reorganized somehow

    static std::vector<std::complex<double> > colStor;
    static UseMatrix* colMat;

    void fuse(); // fuse the compatible children (i.e. same iIndex,jIndex, absorbable oFloor
    void fuseBottomUp(); // fuse tree bottom up
    /// if isLeaf(), but has not block try apply's from  derived class(es)
    void applyDerivedClass(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    void addMulti(const std::string Name, const OperatorDefinition & Definition, const Index *ISub, const Index *JSub,
                  std::complex<double> Multiplier, std::vector<std::complex<double> *> TFac);
    bool absorb(OperatorTree *Other); // absorb Other into present, return false if cannot be absorbed
    static bool skipAxis(const Index* Idx); // true where axis will not be counted as operator-axis

    /// (obsolete)
    void _matrixContractedOLD(UseMatrix &Mat) const;
    /// construct matrix with continuity conditions imposed (columns/rows contracted), exploit block-structure
    void _matrixContracted(UseMatrix &Mat, const Index* IRoot, const Index* JRoot, const std::vector<unsigned int>& ICont, const std::vector<unsigned int>& JCont) const;
    void _blockContracted(UseMatrix &Mat, int & I0, int & J0, const Index* IRoot, const Index* JRoot,
                          const std::vector<unsigned int> &ICont, const std::vector<unsigned int> &JCont) const; ///< contracted block and global position

    void _buildFromMatrix(const Index* Idx, const Index* Jdx, const Eigen::MatrixXcd &Mat);

    /// main internal operator constructor: also holds pointers to time-dependent factors and constant multipliers
    OperatorTree(const std::string Name, const OperatorDefinition & Definition, const Index* IIndex, const Index* JIndex,
                 std::complex<double> Multiplier, std::vector<std::complex<double> *> TFac);

    /// auxiliary internal constructor
    OperatorTree(std::string Name, const Index* Idx, const Index* Jdx, const Eigen::MatrixXcd& Mat);

    OperatorTree(const OperatorTree& Other); // suppress copy constructor

//    // disable Tree::purge for outside use (suppresses warnings about hiding by other purges here)
//    void purge(unsigned int Height) override;

protected:
    //HACK both of these should really be constant
    OperatorFloor * oFloor;
    bool _view;

    void _apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    bool isSymmetric(std::string Kind, double Eps) const override;

public:
    static bool debug;
    OperatorTree():oFloor(nullptr),_view(false){}
    OperatorTree(std::string Name,const Index* IIndex,const Index* JIndex):OperatorAbstract(Name,IIndex,JIndex),oFloor(nullptr),_view(false){}

    /// main constructor
    OperatorTree(const std::string Name /** free name string */, const std::string & Definition /** valid operator definition */,
                 const Index* IIndex /** lhs discretization */, const Index* JIndex /** rhs discretization */);

    /// auxiliary constructor for OperatorDiagsPermuted, OperatorRALL
    OperatorTree(const std::string Name, const std::string Definition, const Index* IIndex, const Index* JIndex, OperatorFloor *OFloor);

    /// auxiliary constructor for InitialState, OperatorExpandIndex, OperatorMapChannelsSurface
    OperatorTree(const OperatorAbstract *A, const Index* IIndex=nullptr, const Index* JIndex=nullptr);

    virtual ~OperatorTree();
    // need copy and asignement...

    /// add Term to present operator and try fuse blocks
    OperatorTree& add(const OperatorTree * Term);

    /// place Block into operator, create branches as needed
    /// upper part of hiearchy is for the column-index
    /// sorting is same as in Index
    OperatorTree& addColumnwise(OperatorTree * Block, bool View=false);

    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const override;
    void applyTranspose(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

    void postProcess(); ///< post-process special floors, parallel, fuse, purge

    void floorInvert(); ///< replaces floor operators with their inverses, ABORTs if not invertible

    /// replace left and/or right index IRep/JRep=0: do not replace
    /// indices must be structurally compatible, but do not need to be equivalent
    void replaceIndex(const Index* IRep, const Index* JRep);

    bool isZero(double Eps=1.e-12, bool Stochastic=true) const override;
    double norm() const override;

    typedef bool(*purgeCriterion)(const OperatorTree* Op);
    using Tree::purge;
    void purge(double Eps=1.e-12); ///< remove branches with isZero(Eps)==true
    void purge(purgeCriterion Crit );

    using Tree::str; // there is also OperatorAbstract::str
    std::string strNode(int Digits=Tree_defaultKind) const override;
    
    void write(std::ofstream &File);
    bool fromCache(std::string Name);
    OperatorTree(std::ifstream &File,std::string Name="", const Index* IIndex=0, const Index* JIndex=0);

    virtual long applyCount() const override;
    
    std::string def() const { return definition.str(); }

    /*
     * Methods to enable rearranging of optrees
     */
    bool nodeEquivalent(const OperatorTree* Other) const override;
    void nodeCopy(const OperatorTree* Other, bool View) override;

    bool isBlockDiagonal() const override ;
    bool isDiagonal() const override ;
    Coefficients  diagonal(bool OnlyForDiagonalOperator) const;

    /// matrix of block-norms of an operator
    UseMatrix matrixBlocks(unsigned int IDepth=INT_MAX,unsigned int JDepth=INT_MAX) const;

    virtual const OperatorFloor* floor() const;
    virtual OperatorFloor*& floor();

    /// reconstruct operator with floor levels up, reset index floors
    void reFloor(size_t IFloor,size_t JFloor);

    /// construct matrix w/o continuity conditions imposed, exploit block-structure
    void matrix(UseMatrix &Mat) const override;

    /// sub-matrix w/o continuity conditions (returns ref to Mat, defaults to complete matrix)
    Eigen::MatrixXcd  & matrix(Eigen::MatrixXcd & Mat, const Index* ISub=0, const Index* JSub=0) const;
    using OperatorAbstract::matrix;

    /// matrix w/o continuity conditions
    Eigen::MatrixXcd matrix() const override;

    /// matrix of norms up to depth MaxDepth or at most floors
    Eigen::MatrixXd matrixNorms(size_t MaxDepth) const;

    /// sparse matrix for operator, near-zeros by Eps are omitted (see BlockView for exact interpretation of Eps)
    Eigen::SparseMatrix<std::complex<double>> matrixSparse(bool Contract=false, double Eps=1.e-12) const override;

    /// construct matrix with continuity conditions imposed (columns/rows contracted), exploit block-structure
    Eigen::MatrixXcd matrixContracted() const override;

    /// create test operators with Dim-dimensional discretization from file (generated if no file is given)
    static OperatorTree* testOp(int Dim, std::string File="");

    /// compare two operators (with tolerance Eps)
    void compare(const OperatorTree * Other, double Eps=0.) const;
    size_t diagnoseSizeOfNode() const override;

    void updateNonLin(double Time, const Coefficients *C);
    void update(double Time, const Coefficients* CurrentVec=nullptr) override;
    const OperatorTree& time(double Time){update(Time);return *this;}
};

#endif // OPERATORTREE_H
