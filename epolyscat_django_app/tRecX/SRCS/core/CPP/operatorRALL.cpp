// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorRALL.h"


#include "qtEigenDense.h"
#include "printOutput.h"

#include "operatorTree.h"
#include "operatorMap.h"
#include "operatordDiagsPermuted.h"
#include "operatorFProduct.h"

#include "blockView.h"
#include "basisVector.h"

#include "indexNew.h"

#include "eigenTools.h"


static OperatorTree* tripletOperator(const BlockView& View, int Triplet){

    OperatorTree* res=new OperatorTree(View.name+tools::str(Triplet),View.block(Triplet,0)->iIndex,View.block(Triplet,0)->jIndex);
    for(int b=0;b<View.blocks(Triplet);b++){
        const OperatorTree* bt=(dynamic_cast<const OperatorTree*>(View.block(Triplet,b)));
        if(not bt)DEVABORT("need OperatorTree type block");
        res->childAdd(bt->detach());
    }
    return res;
}

/// select the MaxU largest singular vectors (with singular value above threshold)
static void addLargeSingular(size_t MaxU, double ThresholdSq,
                             std::vector<std::complex<double>> &SRow,std::vector<Eigen::VectorXcd> &URow,
                             const Eigen::VectorXcd &Val, const Eigen::MatrixXcd &Vec){
    for (int k=0;k<Vec.cols();k++){
        if(std::norm(Val(k))<ThresholdSq)continue;

        size_t ins=SRow.size();
        while(0<ins and std::norm(SRow[ins-1])<std::norm(Val(k)))ins--;
        SRow.insert(SRow.begin()+ins,Val(k));
        URow.insert(URow.begin()+ins,Vec.col(k));
    }

    // truncate, but count similar values only once
    int nSame=0;
    for(size_t i=1;i<SRow.size();i++)nSame+=int(std::norm(SRow[i-1])-std::norm(SRow[i])<std::abs(SRow[0])*1.e-5);
    SRow.resize(std::min(SRow.size(),MaxU+nSame));
    URow.resize(SRow.size());
}

static Eigen::MatrixXcd schmidtOrthonormalize(std::vector<Eigen::VectorXcd> &URow,size_t Truncate){
    double Eps=1.e-12;
    std::vector<int> nonz(URow.size(),0);
    size_t cntn=0;
    for(size_t i=0;i<URow.size() and cntn<Truncate; i++){
        for(size_t j=0;j<i;j++)
            URow[i]-=URow[j]*(URow[j].adjoint()*URow[i]);

        std::complex<double> normI=URow[i].adjoint()*URow[i];
        if(std::abs(normI)<Eps)URow[i].setZero();
        else                  {URow[i]/=sqrt(normI);nonz[i]=1;cntn++;}
    }
    Eigen::MatrixXcd u;
    if(URow.size()>0){
        u=Eigen::MatrixXcd(URow.front().size(),cntn);
        for(size_t i=0,j=0;i<URow.size() and j<cntn;i++){
            if(nonz[i]){u.col(j)=URow[i];j++;}
        }
    }
    return u;
}

/// square root of inverse overlap, weighted for section size:
///     weig/sqrt(diag[Idx])
/// weig = sum[i] |diag[i]| (1-norm of diag)
static void addOvrWeig(std::map<const Index*,Eigen::VectorXcd>&OvrWeig,const Index* Idx){
    if(OvrWeig.count(Idx))return;

    //HACK really an ugly hack...
    const_cast<Index*>(Idx->root())->localOverlapAndInverse(0,0);

    Eigen::MatrixXcd ovI=Idx->localOverlap()->matrix();
    // overlaps must be diagonal
    if(not ovI.isDiagonal())DEVABORT(Sstr+"overlap is not locally diagonal\n"+ovI);
    Eigen::VectorXcd diag=ovI.diagonal();
    double weig=diag.lpNorm<1>();
    OvrWeig[Idx]=Eigen::VectorXcd((weig*diag.cwiseInverse()).cwiseSqrt());
}

OperatorRALL::OperatorRALL(OperatorTree& O, std::vector<size_t> Height, size_t BlockBand, size_t MaxTerms)
    :OperatorTree(O.name+"_rall",O.iIndex,O.jIndex),maxTerms(MaxTerms)
{
    double eps=1.e-10;

    if(not Height.size()){
        // recursion terminated: move all children from input to new
        OperatorTree* op=new OperatorTree(O.name,O.iIndex,O.jIndex);
        O.childrenMove(*op);
        _terms.push_back(_term(op,0,0));
        return;
    }

    // get a block-view at level Height[0] above floor
    BlockView view(&O,O.iIndex->heightAboveFloor()-Height[0],O.jIndex->heightAboveFloor()-Height[0]);


    if(O.idx()!=O.jdx())DEVABORT("must have identical lhs and rhs index pointers");
    if(not O.isComplexSymmetric())DEVABORT(EigenTools::str(O.matrix(),1)+"\nfor now only for complex symmetric operators");

    // move banded part of O to first term (if desired), point to the rest from "residual"
    std::map<const Index*,std::map<const Index*,std::unique_ptr<OperatorTree>>>residual; // matrix blocks of full operator - current approx
    OperatorTree* band(0);
    if(not _terms.size() and BlockBand>0)band=new OperatorTree("band_rall",O.idx(),O.jdx());
    std::map<const Index*,Eigen::VectorXcd>ovrWeig;
    for(int t=0;t<view.triplets();t++){
        OperatorTree* op=tripletOperator(view,t);
        addOvrWeig(ovrWeig,op->iIndex);
        addOvrWeig(ovrWeig,op->jIndex);
        if(band and std::abs(view.tripletI(t)-view.tripletJ(t))<int(2*BlockBand)-1)
            band->addColumnwise(op);
        else
            residual[op->iIndex][op->jIndex].reset(op);
    }
    if(band)_terms.push_back(_term(band,0,0));

    OperatorTree::debug=false;

    while(_terms.size()<maxTerms){
        // accumulate largest singular vectors for present term
        std::map<const Index*,std::map<const Index*,Eigen::VectorXcd>>diags;
        std::map<const Index*,std::vector<Eigen::VectorXcd>>uBlock;
        std::map<const Index*,std::vector<std::complex<double>>> sRow;
        bool allDiagonals=true;
        for(auto &row: residual){
            for(auto& block: row.second){
                if(row.first>block.first)continue; // work with lower triangle

                Eigen::MatrixXcd mat=block.second->matrix();

                // collect diagonals, if all are diagonal
                if(allDiagonals and mat.isDiagonal()){
                    Eigen::VectorXcd d=mat.diagonal();
                    if(not d.isZero()){
                        diags[row.first][block.first]=d;
                        diags[block.first][row.first]=d;
                    }
                }
                else
                    allDiagonals=false;

                // get S M S, assuming S is diagonal
                mat=ovrWeig[block.second->iIndex].asDiagonal()*mat*ovrWeig[block.second->jIndex].asDiagonal();

                // get block singular values (or eigenvalues)
                if(block.second->idx()!=block.second->jdx()){
                    Eigen::BDCSVD<Eigen::MatrixXcd>svd(mat,Eigen::ComputeThinU|Eigen::ComputeThinV);
                    // add to set of largest, possibley removing smaller ones
                    addLargeSingular(mat.rows(),eps,sRow[block.second->iIndex],uBlock[block.second->iIndex],svd.singularValues(),svd.matrixU());
                    addLargeSingular(mat.cols(),eps,sRow[block.second->jIndex],uBlock[block.second->jIndex],svd.singularValues(),svd.matrixV().conjugate());
                }
                else {
                    // on diagonal, ensure symmetric SVD, i.e. standard spectral decomp
                    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
                    ces.compute(mat);
                    addLargeSingular(mat.rows(),eps,sRow[block.second->iIndex],uBlock[block.second->iIndex],ces.eigenvalues(),ces.eigenvectors());
                }

            }
        }

        if(allDiagonals){
            // matrix of diagonal blocks - resort and iterate the scheme
            if(not diags.size())break; // operator is Zero

            OperatorDiagsPermuted *odia=new OperatorDiagsPermuted(O.name+"_diagBlock",O.iIndex,O.jIndex,diags);
            if(odia->isEfficient()){
                // iterate only if scheme is indeed efficient
                if(not odia->symmetrize())DEVABORT("DiagsPermuted is not symmetric");
                OperatorTree *d=odia->diag;
                odia->diag=new OperatorRALL(*d,{0},0,3);
                delete d;
                _terms.push_back(_term(odia,0,0));
                residual.clear();
                break;
            }
            else
                delete odia;
        };
        // all transformations complete

        // Schmidt-orthonormalize and place into Eigen::MatrixXcd
        std::map<const Index*,Eigen::MatrixXcd> uTrans;
        for(auto &u: uBlock){
            if(u.second.size()>0){
                uTrans[u.first]=schmidtOrthonormalize(u.second,u.first->size());
                u.second.clear();
            }
        }

        // no new vectors - terminate terms
        if(uTrans.size()==0)break;

        // build the singular value index
        Index *singI(new Index(*O.iIndex));
        for(auto t: uTrans)singI->nodeAt(t.first->parent()->index())
                ->childReplace(t.first->nSibling(),new Index({new BasisVector(t.second.cols())},{"Vec"}));
        singI->sizeCompute();
        Index *singJ=singI; // for now, symmetric problems only

        // construct new term
        std::map<const Index*,const OperatorFloor*> uFloor,vFloor;
        OperatorTree *uOp=new OperatorTree("U",O.idx(),singI);
        OperatorTree *vOp=new OperatorTree("Vt",singJ,O.jdx());
        for (auto &u: uTrans){
            // construct left (uT) and right (vT) block-wise transformations
            Eigen::MatrixXcd uT=ovrWeig[u.first].cwiseInverse().asDiagonal()*uTrans[u.first];
            uFloor[u.first]=OperatorFloor::factory(std::vector<const Eigen::MatrixXcd*>(1,const_cast<const Eigen::MatrixXcd*>(&uT)),"U");
            const Index* bI=singI->nodeAt(u.first->index());
            uOp->addColumnwise(new OperatorTree("U"+tools::str(u.first->levelRank()),"n/a",O.idx()->nodeAt(u.first->index()),bI,
                                                const_cast<OperatorFloor*>(uFloor[u.first])));
            Eigen::MatrixXcd vT=uT.transpose();
            vFloor[u.first]=OperatorFloor::factory(std::vector<const Eigen::MatrixXcd*>(1,&vT),"Vt");
            vOp->addColumnwise(new OperatorTree("Vt"+tools::str(u.first->levelRank()),"n/a",bI,O.jdx()->nodeAt(u.first->index()),
                                                const_cast<OperatorFloor*>(vFloor[u.first])));
        }

        std::string n=tools::str(_terms.size());
        OperatorTree *qDiags= new OperatorTree("Si"+n,singI,singJ); // need operator tree access
        for(auto & row: residual){
            Str wNameI("Q_","");
            for(auto & block: row.second){
                const Index* blockI=block.second->iIndex;
                const Index* blockJ=block.second->jIndex;
                Index* blockU=singI->nodeAt(blockI->index());
                Index* blockV=singJ->nodeAt(blockJ->index());

                if(not uTrans.count(blockI) or not uTrans.count(blockJ))continue; // no bases on this block
                if(not blockU->size() or not blockV->size())continue;

                Eigen::MatrixXcd mat(block.second->matrix());
                Eigen::MatrixXcd &uI=uTrans[blockI],&vJ=uTrans[blockJ];

                // get diagonal values in transformed representation
                Eigen::MatrixXcd diag=Eigen::MatrixXcd::Zero(uI.cols(),vJ.cols());
                for(int i=0;i<std::min(diag.rows(),diag.cols());i++)diag(i,i)=
                        ((uI.col(i).conjugate().cwiseProduct(ovrWeig[blockI])).transpose()
                         *mat*(vJ.col(i).conjugate().cwiseProduct(ovrWeig[blockJ])))(0);

                // insert diagonal floor
                OperatorFloor* wFloor=OperatorFloor::factory(std::vector<const Eigen::MatrixXcd*>(1,&diag),"uTrans");
                qDiags->addColumnwise(new OperatorTree(wNameI+block.first,"n/a",blockU,blockV,wFloor));

                // update residual: subtract uT[i,a] Q[a,a'] vT[i,a] from resdiual operator (views!)
                if(uFloor.count(row.first) and vFloor.count(block.first)){
                    OperatorFProduct* rf=new OperatorFProduct(-1.,{uFloor[row.first],wFloor,vFloor[block.first]},true);
                    residual[row.first][block.first]
                            ->childAdd(new OperatorTree(wNameI+block.first->levelRank(),"n/a",blockI,blockJ,rf));
                }

            }
        }
        if(qDiags->childSize()){
            // recursively analyze w, if Height is given for next recursive level
            _terms.push_back(_term(new OperatorRALL(*qDiags,{Height.begin()+1,Height.end()},MaxTerms),uOp,vOp));
            if(_terms.size()==MaxTerms)PrintOutput::warning(Sstr+"exhausted maximal number of"+MaxTerms
                                                            +"terms in OperatorRALL, approximation may not be accurate to "+eps);
        }
        else {
            delete uOp;delete vOp;delete qDiags;
        }
    }
}

void OperatorRALL::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    Y.scale(B);
    for(const _term &t: _terms){
        if(t.u and t.v){
            t.v->apply(A,Vec,0.,*t.tmpR);
            t.si->apply(1.,*t.tmpR,0.,*t.tmpL);
            t.u->apply(1.,*t.tmpL,1.,Y);
        }
        else if (not t.u and not t.v){
            t.si->apply(A,Vec,1.,Y);
        }
        else
            if(not(t.u and t.v))DEVABORT("must have both or no transformation");
    }
}

