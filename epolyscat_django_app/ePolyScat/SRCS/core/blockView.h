// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BLOCKVIEW_H
#define BLOCKVIEW_H

#include <vector>
#include <complex>

#include "operatorAbstract.h"
#include "index.h"
#include "operatorAbstract.h"

#include "operatorTree.h"
#include "qtEigenSparse.h"


///@brief view operator as a sparse matrix of Operators
class BlockView: public OperatorAbstract
{
private:
    struct Triplet{
        int i,j;
        Triplet(){}
        Triplet(int I, int J):i(I),j(J){}
        std::vector<const OperatorAbstract*> block;
        void add(const OperatorAbstract *Block);
        int blockRow() const {return i;}
        int blockCol() const {return j;}
    };
    std::vector<Triplet> _trip;

    void construct(const OperatorTree* Op, const int IDepth, const int JDepth,
                   std::map<const Index*,int> &ICache,
                   std::vector<std::vector<std::pair<int,int>>>& Loc
                   );
    void renumber(); // introduce consistent numbering for blocks

    int superDiagonals() const; ///< number of superdiagonals as for banded storage
    int subDiagonals() const; ///< number of subdiagonals as for banded storage

    /// returns one BlockView for each multiplicative factor (for now: only if floor blocks)
    void splitByFactors(std::vector<BlockView> & Part, std::vector<std::complex<double> *> &Factor);


//    void update(double Time, const Coefficients* CurrentVec){} /// dummy (no time-updates for BlockView)
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

public:
    static const int depthFloor=INT_MAX;
    static const int depthContinuity=INT_MAX-1;

    virtual ~BlockView(){}
    BlockView(){}
    BlockView(const OperatorTree *Op, int IDepth=depthFloor, int JDepth=depthFloor);
    BlockView(const OperatorAbstract *Op);

    int blockRows() const;
    int blockCols() const;
    int tripletI(int Triplet){return _trip[Triplet].i;};
    int tripletJ(int Triplet){return _trip[Triplet].j;};

    std::vector<std::complex<double> > bandedStorage() const; ///< data in banded storage (LAPACK style)
    /// matrix in sparse storage: mat(i,j) <= Eps sqrt( max|row(i)|* max|col(j)| ) are omitted
    void sparseMatrix(Eigen::SparseMatrix<std::complex<double> > & Mat, bool Contract=false /** */, double Eps=0.) const;

    std::string normsMatrix(int Digits=0);
    std::string str(std::string Info="") const;
    void printRowsCols() const; ///< table of row and block indices

    int triplets() const{return _trip.size();} ///< number of block-triplets
    int blocks(int Triplett) const {return _trip[Triplett].block.size();} ///< number of blocks on triplets
    const OperatorAbstract* block(int Triplett, int Block) const {return _trip[Triplett].block[Block];} ///< Block's block on triplet

    bool isSymmetric(std::string Kind, double Eps) const;
};

#endif // BLOCKVIEW_H
