// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "blockView.h"

#include <utility>

#include "timer.h"
#include "str.h"
#include "index.h"
#include "operatorFloor.h"
#include "parameters.h"
#include "printOutput.h"

#include "log.h"
#include "treeIterator.h"

#include "mpiWrapper.h"
#include "parallelOperator.h"

using namespace std;

void BlockView::apply(complex<double> A, const Coefficients &Vec,  complex<double> B, Coefficients &Y) const
{
    ABORT("not implemented");
}

int cacheIndexBlocks(map<const Index*, std::pair<int, std::size_t>> Cache, int Depth, const Index* Idx){
    // identify block indices: either hasFloor() ABOVE Depth or block->depth()==Depth
    int nBlocks(0);
    for(const Index* ix=Idx;ix!=0;ix=ix->nodeNext()){
        int d=ix->depth();
        if(d>=Depth and (ix->hasFloor() or d==Depth)){
            Cache[ix]=std::make_pair(nBlocks,ix->depth());
            nBlocks++;
        }
    }
    return nBlocks;
}

void BlockView::Triplet::add(const OperatorAbstract *Block){
    block.push_back(Block);
}
BlockView::BlockView(const OperatorAbstract * Op)
    :OperatorAbstract("view("+Op->name+")",Op->iIndex,Op->jIndex){
    _trip.push_back(Triplet(0,0));
    _trip.back().block.push_back(Op);
}

int getIndexBlock(const Index* Idx,std::map<const Index*,int> &Cache ){
    // number blocks by first request
    // note: this may create a funny numbering of the blocks
    auto ic=Cache.find(Idx);
    int iBlock;
    if(ic==Cache.end()){
        iBlock=Cache.size();
        Cache[Idx]=iBlock;
    }
    else
        iBlock=ic->second;
    return iBlock;
}

void BlockView::construct(const OperatorTree* Op,const int IDepth, const int JDepth,
                          std::map<const Index*,int> &ICache,
                          std::vector<vector<pair<int,int>>>& OpBlocks
                          ){


    if(Op->floor() or
            (       Op->iIndex->depth()==IDepth and Op->child(0)->iIndex->depth()!=size_t(IDepth)
                    and Op->jIndex->depth()==JDepth and Op->child(0)->jIndex->depth()!=size_t(JDepth))){
        // add block to sparse (stored as Tripletts):
        // if iBlock,jBlock pair did not occur previously
        //    create a new triplet of iBlock,jBlock,list-of-terms
        // append block to lists-of-terms

        // OpBlock: points from index-pair to triplets as OpBlocks[i][j].first=j

        // get block number from cache or create new
        // for consistent block numbering we should create the complete cache
        // in a first run, sort index, and re-run
        size_t iBlock=getIndexBlock(Op->iIndex,ICache);
        int jBlock=getIndexBlock(Op->jIndex,ICache);

        bool newBlock=true;
        if(OpBlocks.size()<=iBlock)OpBlocks.resize(iBlock+1);
        for(auto& v: OpBlocks[iBlock]){
            if(v.first==jBlock){ // search through block index pairs
                _trip[v.second].block.push_back(Op);
                newBlock = false;
                break;
            }
        }

        if(newBlock){
            _trip.push_back(Triplet(iBlock,jBlock)); // append new block index pair
            OpBlocks[iBlock].push_back({jBlock,_trip.size()-1}); // set pointer to new pair
            _trip.back().block.push_back(Op); // given iBlock,jBlock
        }
    }
    else {
        // Op not accepted, proceed down from node
        for (size_t k=0;k<Op->childSize();k++)
            construct(Op->child(k),IDepth,JDepth,ICache,OpBlocks);
    }
}

void nodeNumbers(std::map<const Index*,int>&Cache,const Index* Beg, size_t Depth){
    for(const Index* ix=Beg;ix!=0;ix=ix->nodeNext(Beg)){
        if((ix->depth()==Depth or ix->hasFloor()) and Cache.count(ix)==0)Cache[ix]=Cache.size();
    }
}

BlockView::BlockView(const OperatorTree * Op,  int IDepth, int JDepth)
    :OperatorAbstract("view("+Op->name+")",Op->iIndex,Op->jIndex)
{
    std::map<const Index*,int> iCache;
    nodeNumbers(iCache,Op->jIndex,JDepth);
    nodeNumbers(iCache,Op->iIndex,IDepth);

    std::vector<vector<pair<int,int>>> opBlocks;
    construct(Op,IDepth,JDepth,iCache,opBlocks);

    // remove empty triplets
    for(int k=_trip.size()-1;k>=0;k--)
        if(!_trip[k].block.size())_trip.erase(_trip.begin()+k);

    // re-number: count only blocks that actually appear
    renumber();
    int nTrips=_trip.size();
    if(ParallelOperator::getHost(Op,false)==ParallelOperator::distributed)
        MPIwrapper::AllreduceMAX(&nTrips,1);
    if(nTrips!=int(_trip.size()))DEVABORT(Sstr+"trips differ"+Op->name+nTrips+_trip.size());
    for(auto t: _trip)
        for(auto b: t.block)
            if(b->iIndex!=t.block.front()->iIndex or b->jIndex!=t.block.front()->jIndex){
                for(auto c:t.block)Sstr+c->iIndex->strNode()+c->jIndex->strNode()+Sendl;
                DEVABORT("blocks incorrectly assembled");
            }
}

void BlockView::renumber(){
    {
        std::set<int>in; // set is guaranteed to be sorted in increasing order
        for(auto t: _trip)in.insert(t.i);
        std::map<int,int> sort;
        for(auto i: in)sort[i]=sort.size();
        for(auto &t: _trip)t.i=sort[t.i];
    }
    {
        std::set<int>in; // set is guaranteed to be sorted in increasing order
        for(auto t: _trip)in.insert(t.j);
        std::map<int,int> sort;
        for(auto i: in)sort[i]=sort.size();
        for(auto &t: _trip)t.j=sort[t.j];

    }
}


string BlockView::normsMatrix(int Digits){
    UseMatrix nrm=UseMatrix::Zero(blockRows(),blockCols());
    for(size_t k=0;k<_trip.size();k++){
        complex<double> a=0.;
        for (size_t l=0;l<_trip[k].block.size();l++){
            a+=_trip[k].block[l]->norm();
        }
        nrm(_trip[k].blockRow(),_trip[k].blockCol())=a;
    }
    return nrm.str("",Digits);
}

int BlockView::subDiagonals() const {
    int d=0;
    for(int i=0;i<_trip.size();i++){
        const OperatorAbstract *b=_trip[i].block[0];
        d= max(d,int(b->jIndex->posIndex()+b->jIndex->sizeStored())-1-int(b->iIndex->posIndex()));
    }
    return d;
}

int BlockView::superDiagonals() const {
    int d=0;
    for(int i=0;i<_trip.size();i++){
        const OperatorAbstract *b=_trip[i].block[0];
        d= max(d,int(b->iIndex->posIndex()+b->iIndex->sizeStored())-1-int(b->jIndex->posIndex()));
    }
    return d;
}

vector< complex<double> > BlockView::bandedStorage() const {
    int col=jIndex->sizeStored(),sup=superDiagonals();
    int ld=subDiagonals()+sup+1;
    vector< complex<double> > stor(ld*col,0.);

    for(size_t k=0;k<_trip.size();k++){
        int i0=_trip[k].block[0]->iIndex->posIndex(),m=_trip[k].block[0]->iIndex->sizeStored();
        int j0=_trip[k].block[0]->jIndex->posIndex(),n=_trip[k].block[0]->iIndex->sizeStored();
        vector< complex<double> > mat;
        for(size_t l=0;l<_trip[k].block.size();l++){
            _trip[k].block[l]->OperatorAbstract::matrix(mat);
            for(int j=0,ji=0;j<n;j++)
                for(int i=0;i<m;i++,ji++)
                    stor[sup+(i0+i-j0+j) + (j0+j)*ld]+=mat[ji];
        }
    }
    return stor;
}


void BlockView::sparseMatrix(Eigen::SparseMatrix< complex<double> > &Mat,bool Contract,double Eps) const {

    vector<Eigen::Triplet< complex<double> > > list;

    int idim,jdim;
    vector<unsigned int> iglob,jglob;
    vector<double>inrm,jnrm;
    if(Contract){
        idx()->multiplicities(iglob,inrm);
        jdx()->multiplicities(jglob,jnrm);
    }
    else {
        for(size_t k=0;k<iIndex->sizeStored();k++)iglob.push_back(k);
        for(size_t k=0;k<jIndex->sizeStored();k++)jglob.push_back(k);
        idim=iglob.size();
        jdim=iglob.size();
    }

    // convert to triples wrt to global index, normalize where multiplicities
    idim=0;
    jdim=0;
    Parameters::updateToOne(); // set all time-dependent parameters =1
    for(size_t k=0;k<_trip.size();k++){
        int i0=_trip[k].block[0]->iIndex->posIndex(iIndex),m=_trip[k].block[0]->iIndex->sizeStored();
        int j0=_trip[k].block[0]->jIndex->posIndex(jIndex),n=_trip[k].block[0]->jIndex->sizeStored();
        vector< complex<double> > mat;
        for(size_t l=0;l<_trip[k].block.size();l++){
            _trip[k].block[l]->matrix(mat);

            std::vector<double> rMax,cMax;
            for(int i=0;i<m;i++)rMax.push_back(Eigen::Map<Eigen::MatrixXcd>(mat.data(),m,n).row(i).lpNorm<Eigen::Infinity>());
            for(int j=0;j<n;j++)cMax.push_back(Eigen::Map<Eigen::MatrixXcd>(mat.data(),m,n).col(j).lpNorm<Eigen::Infinity>());

            for(int j=j0,ij=0;j<j0+n;j++){
                double epsSqJ=Eps*Eps*cMax[j-j0];
                for(int i=i0;i<i0+m;i++,ij++){
                    complex<double> matij=mat[ij];
                    if(inrm.size()>0)matij*=inrm[i];
                    if(jnrm.size()>0)matij*=jnrm[j];
                    idim=max(idim,(int)iglob[i]);
                    jdim=max(jdim,(int)jglob[j]);
                    if(std::norm(matij)>epsSqJ*rMax[i-i0])list.push_back(Eigen::Triplet< complex<double> >(iglob[i],jglob[j],matij));
                }
            }
        }
    }
    Parameters::restoreToTime();
    idim++;
    jdim++;

    Mat.resize(idim,jdim);
    Mat.setFromTriplets(list.begin(),list.end(),[] (const  complex<double> & a,const  complex<double> &b) { return a+b; });
}

void BlockView::splitByFactors(std::vector<BlockView> &Part, std::vector<std::complex<double>* > &Factor){

    for(int k=0;k<_trip.size();k++){
        for(int l=0;l<_trip[k].block.size();l++){
            const OperatorTree* op=dynamic_cast<const OperatorTree*>(_trip[k].block[l]);
            if(op==0 or op->floor()==0)ABORT(Str("only implemented for floor blocks")+op);
            int iPart=std::find(Factor.begin(),Factor.end(),op->floor()->factor())-Factor.begin();
            if(iPart>=Factor.size()){
                // create index structure, but remove all blocks
                Part.push_back(*this);
                for(int m=0;m<_trip.size();m++)Part.back()._trip[m].block.clear();
                Factor.push_back(op->floor()->factor());
            }
            Part[iPart]._trip[k].add(op); // add to matching part
        }
    }

    // remove all triplets that have remained empty
    for(size_t k=0;k<Part.size();k++){
        for(int l=Part[k]._trip.size()-1;l>=0;l--){
            if(Part[k]._trip[l].block.size()==0)
                Part[k]._trip.erase(Part[k]._trip.begin()+l);
        }
    }
}

int BlockView::blockRows() const {
    int n=0;
    for(size_t k=0;k<_trip.size();k++)n= max(n,_trip[k].blockRow());
    return n+1;
}

int BlockView::blockCols() const {
    int n=0;
    for(size_t k=0;k<_trip.size();k++)n= max(n,_trip[k].blockCol());
    return n+1;
}

string BlockView::str(string Info) const {
    Str s;
    for(size_t k=0;k<_trip.size();k++)
    {
        s+=Sstr+tools::str(_trip[k].i)+" "+tools::str(_trip[k].j)+_trip[k].block.front()->iIndex->strNode()+_trip[k].block.front()->jIndex->strNode();
        for(size_t l=0;l<_trip[k].block.size();l++)s+=_trip[k].block[l]->str()+" \n ";
    }
    return std::move(s);
}

void BlockView::printRowsCols() const {
    std::vector<std::string> rows;
    std::vector<std::string> cols;
    for(auto t: _trip){
        if((int) rows.size()<=t.i)rows.resize(t.i+1);
        if((int) cols.size()<=t.j)cols.resize(t.j+1);
        rows[t.i]=t.block[0]->idx()->strNode();
        cols[t.j]=t.block[0]->jdx()->strNode();
    }
    PrintOutput::newRow();
    PrintOutput::rowItem("#");
    PrintOutput::rowItem("rowIdx");
    PrintOutput::rowItem("colIdx");
    for(size_t k=0;k< std::max(rows.size(),cols.size());k++){
        PrintOutput::newRow();
        PrintOutput::rowItem(int(k));
        PrintOutput::rowItem((k<rows.size()?rows[k]:"-"));
        PrintOutput::rowItem((k<cols.size()?cols[k]:"-"));
    }
    PrintOutput::end();
}

bool BlockView::isSymmetric(std::string Kind, double Eps) const{
    for(auto t: _trip){
        if(t.i<= t.j){
            // get matching triplet
            for(auto tt: _trip){
                if(tt.i==t.j and tt.j==t.i){
                    if(tt.block.size()!=t.block.size())return false;
                    if(tt.block.size()==0)continue;

                    if(     (tt.block.front()->iIndex!=t.block.front()->jIndex or
                             tt.block.front()->jIndex!=t.block.front()->iIndex)
                            and
                            (not tt.block.front()->iIndex->treeEquivalent(t.block.front()->jIndex) or
                             not tt.block.front()->jIndex->treeEquivalent(t.block.front()->iIndex) )){
                        Sstr+t.i+t.j+tt.i+tt.j
                                +"\n"+ t.block.front()->iIndex+ t.block.front()->jIndex
                                +"\n"+tt.block.front()->jIndex+tt.block.front()->iIndex
                                +"\n"+ t.block.front()->iIndex->strNode()+ t.block.front()->jIndex->strNode()
                                +"\n"+tt.block.front()->iIndex->strNode()+tt.block.front()->jIndex->strNode()+Sendl;
                        DEVABORT(Index::failureCompatible+"\nnot symmetric - incorrect triplet indices");
                    }

                    Eigen::MatrixXcd bm=Eigen::MatrixXcd::Zero(t.block.front()->iIndex->size(),t.block.front()->jIndex->size());
                    for(auto b: t.block)bm+=b->matrix(); // add all blocks in triplet
                    double nrm=bm.lpNorm<Eigen::Infinity>();
                    if     (Kind=="ComplexSymmetric")for(auto bt: tt.block)bm-=bt->matrix().transpose();
                    else if(Kind=="SelfAdjoint"     )for(auto bt: tt.block)bm-=bt->matrix().adjoint();
                    else DEVABORT("isSymmetric needs Kind=ComplexSymmetric or SelfAdjoint, got: "+Kind);

                    if(not bm.isZero(std::max(nrm,1.)*Eps))return false;
                }
            }
        }
    }
    return true;
}



