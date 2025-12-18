// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef PATCHRADIAL_H
#define PATCHRADIAL_H

#include "qtEigenDense.h"
#include <vector>
#include "densityMatrix1.h"
#include "index.h"
#include "abort.h"
#include "basisIntegrable.h"

class BasisIntegrable;
class OperatorTree;
class OperatorFloor;

///@brief collects operator leaf's of a given Floor type by the Floors' radial bases
///
/// sort into a hierarchy of LeafM[ma,mb].LeafL[la,lb].LeafC[chanA,chanB]
template<class Floor>
class PatchRadial{
    const BasisIntegrable *_basA,*_basB;

    /// get channel number and M,L quantum numbers for given Idx
    static void getChanML(const Index* Idx,int & Chan, int &M, int &L){
        M=INT_MAX;
        L=INT_MAX;
        Chan=INT_MAX;
        while(Idx->parent()!=0){
            if(Idx->parent()->axisName()=="Eta")L=Idx->parent()->basis()->physical(Idx->nSibling());
            if(Idx->parent()->axisName()=="Phi")M=Idx->parent()->basis()->physical(Idx->nSibling());
            if(Idx->parent()->axisName()=="Channel")Chan=Idx->nSibling();
            else if(Idx->parent()->axisName().find("Orbital&")==0)Chan=0;
            if(Chan!=INT_MAX and M!=INT_MAX and L!=INT_MAX){
                if(std::abs(M)>L)DEVABORT(Sstr+"BAD: c,m > l"+Chan+M+L);
                return;
            }
            Idx=Idx->parent();
        }
        ABORT("need axes (Channel or Orbital&..), Eta, and Phi, found hierarchy: "+Idx->root()->hierarchy());
    }

    /// true if Leaf is to be selected into current patch (use present, if no basis pair given)
    bool select(const OperatorTree &Leaf, const BasisIntegrable *&BasA, const BasisIntegrable *&BasB){
        const Floor* f=dynamic_cast<const Floor*>(Leaf.floor());
        // if the Floor has already been set up, do not select it
        if(not f or f->hasBeenSetUp())return false;

        if(BasA==0){
            // when first entering, patch has not been decided upon: BasA=0,BasB=0
            // set BasA and BasB, use patch on first floor that that has not been set up yet
            BasA=Leaf.iIndex->basis()->integrable();
            BasB=Leaf.jIndex->basis()->integrable();
            if(BasA==0 or BasB==0)DEVABORT("must have BasisIntegrable on operator leafs");
            return true;
        }
        else{
            // select Leaf, if radial patch matches patch of BasA x BasB
            return  BasA->lowBound()==Leaf.iIndex->basis()->integrable()->lowBound() and
                    BasB->lowBound()==Leaf.jIndex->basis()->integrable()->lowBound();
        }
    }
public:

    /// patch of all Operator leafs that match first undefined
    const BasisIntegrable* aBas() const {return _basA;} ///< left radial basis of patch
    const BasisIntegrable* bBas() const {return _basB;} ///< right radial basis of patch

    struct LeafC{
        /// matrix elements for a given radial patch and la,ma,lb,mb and channel pair CaCb (channel numbering as in class DensityMatrix1)
        LeafC(OperatorTree* Op, int CaCb, std::string Kind):op(Op),cc(CaCb)
        {
            mat=&dynamic_cast<Floor*>(op->floor())->mat();
            if(Kind=="full")         *mat=Eigen::MatrixXcd::Zero(op->iIndex->size(),op->jIndex->size());
            else if(Kind=="diagonal")*mat=Eigen::MatrixXcd::Zero(op->iIndex->size(),1);
            else ABORT("temporary matrix either Kind=\"full\" or -\"diagonal\"");
        }
        int cc;
        OperatorTree* op; ///< pointer to operator leaf
        Eigen::MatrixXcd *mat; ///< pointer to _mat temporary storage for accumulating matrix
    };

    struct LeafL{
        /// all channel leafs for given ma,mb,la,lb
        LeafL(int La,int Lb):la(La),lb(Lb){}
        int la,lb;
        std::vector<LeafC> leafs;
    };

    struct LeafM{
        /// all channel leafs for given ma,mb
        LeafM(int Ma,int Mb):ma(Ma),mb(Mb){}
        int laMax() const {int ll=0; for(LeafL l: leafs)ll=std::max(ll,l.la); return ll;}
        int lbMax() const {int ll=0; for(LeafL l: leafs)ll=std::max(ll,l.lb); return ll;}

        int ma,mb;
        std::vector<LeafL> leafs;
    };
    std::vector<LeafM> leafs; // all leafs, collected by ma,mb

    /// RadialPatch in Op, using first leaf of Op that not hasBeenSetUp()
    PatchRadial(OperatorTree* Op, const DensityMatrix1 & Rho, std::string Kind)
        :_basA(0),_basB(0){
        // note: rather that running through a constructor, we could do this with some
        // bool PatchRadial::next(Op,Rho,Kind) function and a while loop
        // argument "Kind" should go and be replaced by a function
        // matrixShape<OperatorHartree>()...diagonal, matrixShape<OperatorFloorXC>()...full

        // int cntMll=0; //PARALL - count mll-contributions, distribute across processes
        for(OperatorTree * oLeaf=Op->firstLeaf();oLeaf!=0;oLeaf=oLeaf->nextLeaf()){
            // probably just test here on _basA,_basB rather than hiding this in select(...)
            if(select(*oLeaf,_basA,_basB)){
                // sort patch into list
                int ma,la,mb,lb,chanA,chanB;
                getChanML(oLeaf->iIndex,chanA,ma,la);
                getChanML(oLeaf->jIndex,chanB,mb,lb);

                int m,l;
                for (m=0;m<leafs.size() and (leafs[m].ma!=ma or leafs[m].mb!=mb);m++);
                if(m==leafs.size())leafs.push_back(LeafM(ma,mb));

                for (l=0;l<leafs[m].leafs.size() and (leafs[m].leafs[l].la!=la or leafs[m].leafs[l].lb!=lb);l++);
                if(l==leafs[m].leafs.size())leafs[m].leafs.push_back(LeafL(la,lb));
                leafs[m].leafs[l].leafs.push_back(LeafC(oLeaf,Rho.cacb(chanA,chanB),Kind));


                //PARALL - assign
                // if(MPIwrapper::Rank()==(cntMll++)%MPIwrapper::Size())
                //     leafs[m].leafs[l].leafs.push_back(LeafC(oLeaf,Rho.cacb(chanA,chanB),Kind));
            }
        }
    }
    int laMax(){int lmax=0; for(LeafM l: leafs)lmax=std::max(lmax,l.laMax());return lmax;} ///<largest la in PatchRadial
    int lbMax(){int lmax=0; for(LeafM l: leafs)lmax=std::max(lmax,l.lbMax());return lmax;} ///< larges lb in PatchRadial
    const LeafC & firstC() const {return leafs[0].leafs[0].leafs[0];}
};


#endif // PATCHRADIAL_H
