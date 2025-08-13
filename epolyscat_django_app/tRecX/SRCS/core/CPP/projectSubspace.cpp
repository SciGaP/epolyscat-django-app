// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "projectSubspace.h"

#include "qtEigenDense.h"
#include "eigenTools.h"
#include "coefficients.h"
#include "index.h"
#include "operatorFloor.h"
#include "printOutput.h"
#include "parallelOperator.h"
#include "parallel.h"
#include "inverse.h"
#include "tree.h"
#include "debugInfo.h"

// duplicated in BasisOrbital!
static Eigen::MatrixXcd allCoefs( std::vector<const Coefficients*> ICoefs){
    if(ICoefs.size()==0)DEVABORT("emtpy vector of Coefficients");
    Eigen::MatrixXcd allMat(ICoefs[0]->size(),ICoefs.size());
    for(size_t k=0;k<ICoefs.size();k++){
        if(ICoefs[k]->orderedData())
            allMat.col(k)=Eigen::Map<Eigen::MatrixXcd>(ICoefs[k]->orderedData(),ICoefs[k]->size(),1);
        else
            allMat.col(k)=Eigen::MatrixXcd::Zero(ICoefs[k]->size(),1);
    }
    return allMat;
}

void buildMap(OperatorTree* Map, std::vector<const Coefficients*> C){

    if((Map->iIndex->hasFloor() or Map->jIndex->hasFloor())
            and std::all_of(C.begin(),C.end(),[](const Coefficients* c){return c->isZero();}))
        return;

    if(Map->iIndex->hasFloor() and Map->jIndex->hasFloor()){
        Eigen::MatrixXcd cMat=C[0]->idx()==Map->iIndex ?
                    allCoefs(C) :
                    allCoefs(C).transpose();
        Map->floor()=OperatorFloor::factory(std::vector<const Eigen::MatrixXcd*>(1,&cMat),
                                            "map:"+Map->iIndex->hash()+"|"+Map->jIndex->hash());
    }
    std::vector<const Coefficients*>c(C);

    // column indices
    if(not Map->jIndex->hasFloor())
        for(size_t n=0;n<Map->jIndex->childSize();n++){
            Map->childAdd(new OperatorTree("map",Map->iIndex,Map->jIndex->child(n)));
            if(C[0]->idx()==Map->jIndex)for(size_t l=0;l<c.size();l++)c[l]=C[l]->child(n); // descend in C, if C are columns
            buildMap(Map->childBack(),c);
        }

    // row indices
    if(not Map->iIndex->hasFloor())
        for(size_t n=0;n<Map->iIndex->childSize();n++){
            Map->childAdd(new OperatorTree("map",Map->iIndex->child(n),Map->jIndex));
            if(C[0]->idx()==Map->iIndex)for(size_t l=0;l<c.size();l++)c[l]=C[l]->child(n); // descend in C, if C are rows
            buildMap(Map->childBack(),c);
        }
    if(Map->root()==Map)Map->purge();
}

ProjectSubspace::ProjectSubspace(std::vector<const Coefficients*> Vectors, std::vector<const Coefficients*> Duals, size_t BeginOrthonormal){
    _construct(Vectors,Duals,BeginOrthonormal);
    definition=OperatorDefinition("(from vectors)");
    PrintOutput::message(Str("Set up projector, dim=")+Vectors.size());
}

ProjectSubspace::ProjectSubspace(std::vector<Coefficients *> Vectors, std::vector<Coefficients *> Duals, size_t BeginOrthonormal){
    if(Vectors.size()>0)EigenSolverAbstract::normalize(Duals,Vectors);
    std::vector<const Coefficients*> cV,cD;
    for(Coefficients* c: Vectors)cV.push_back(c);
    for(Coefficients* c: Duals)  cD.push_back(c);
    _construct(cV,cD,BeginOrthonormal);
    definition=OperatorDefinition("(from vectors)");
}


std::vector<std::vector<std::vector<const Coefficients*>>> groupByNonzeros(std::vector<const Coefficients*> Vecs,std::vector<const Coefficients*> Duals,std::vector<int> &Sorting){
    // first nonzero and end of non-zero floors in Coefficient
    struct vStruc {
        int nOrig,beg,end;
        const Coefficients *c,*d;
        vStruc(int NOrig, const Coefficients* C,const Coefficients* D):nOrig(NOrig),c(C),d(D){
            beg=0;
            const Coefficients* f=c->firstLeaf();
            for(;f!=0 and f->isZero();f=f->nextLeaf())beg++;
            for(int k=beg;f!=0;k++,f=f->nextLeaf())
                if(not f->isZero())end=k+1;
        }
    };
    std::vector<vStruc> vS;
    for(size_t k=0;k<Vecs.size();k++)vS.push_back(vStruc(k,Vecs[k],Duals[k]));

    std::stable_sort(vS.begin(),vS.end(),[](const vStruc & A,const vStruc & B){return A.beg<B.beg;});

    std::vector<std::vector<std::vector<const Coefficients*> > > gVecs(2);
    Sorting.clear();
    Sorting.push_back(vS[0].nOrig);
    gVecs[0]=std::vector<std::vector<const Coefficients*> >(1,std::vector<const Coefficients*>(1,vS[0].c));
    gVecs[1]=std::vector<std::vector<const Coefficients*> >(1,std::vector<const Coefficients*>(1,vS[0].d));
    for(size_t k=1;k<vS.size();k++){
        if(vS[k-1].beg<vS[k].beg){
            gVecs[0].push_back(std::vector<const Coefficients*>(0));
            gVecs[1].push_back(std::vector<const Coefficients*>(0));
        }
        Sorting.push_back(vS[k].nOrig);
        gVecs[0].back().push_back(vS[k].c);
        gVecs[1].back().push_back(vS[k].d);
    }
    return gVecs;
}

OperatorTree* transfer(OperatorTree * Op){
    OperatorTree* op=new OperatorTree();
    op->nodeCopy(Op,false);
    Op->childrenMove(*op);
    return op;
}
ProjectSubspace::ProjectSubspace(EigenSolver & Slv){

    if(Slv.isLeaf()){
        // _construct needs vector<const Coefficients*>
        std::vector<const Coefficients*> v,d;
        for(auto c: Slv.rightVectors())v.push_back(c);
        for(auto c: Slv.dualVectors())d.push_back(c);
        _construct(v,d,0);
        return;
    }

    iIndex=Slv.oper()->idx();
    jIndex=Slv.oper()->jdx();

    _subspaceIndex=new Index();
    _subspaceIndex->setAxisName(Str("Subspace","")+Slv.depth());
    _subspaceIndex->setBasis(BasisAbstract::factory(Sstr+"Vector: "+Slv.childSize()));
    _mapTo.reset(new OperatorTree("mapToFullFromContracted",iIndex,_subspaceIndex));
    _mapFrom.reset(new OperatorTree("mapFromFullToContracted",_subspaceIndex,jIndex));

    for(size_t k=0;k<Slv.childSize();k++){
        if(Slv.child(k)->rightVectors().size()==0)continue;

        ProjectSubspace ps(*Slv.child(k));

        // need to "transfer", as _mapFrom will become deleted when ps goes out of scope
        _mapFrom->childAdd(transfer(ps._mapFrom.get()));
        _mapTo->childAdd(transfer(ps._mapTo.get()));
        _subspaceIndex->childAdd(ps._subspaceIndex);
        _sorting.insert(_sorting.end(),ps.sorting().begin(),ps.sorting().end());
    }

    // need to distignuisch hierarchy levels, else it is mistaken for continuity levels
    _subspaceIndex->sizeCompute();
    _subspaceIndex->unsetFloor();
    _mapFrom->purge();
    _mapTo->purge();

    if(Slv.parent()==0){

        // for now, we consider projections as global (can be improved)
        _subspaceC.reset(new Coefficients(_subspaceIndex));


        ParallelOperator::setDistribution(_mapTo.get());
        ParallelOperator::setDistribution(_mapFrom.get());
        PrintOutput::message(Str("Set up projector, dim=")+_subspaceIndex->size());

        // make sure loads are sync'd
        ParallelOperator::sync(_mapFrom.get());
        ParallelOperator::sync(_mapTo.get());
    }
    definition=OperatorDefinition(Slv.selection()+" of "+Slv.oper()->def());
}

static void checkDuals(std::vector<const Coefficients *> Vectors, std::vector<const Coefficients *> Duals,double Epsilon){
    PrintOutput::outputLevel("restore");

    if(Vectors.size()!=Duals.size())DEVABORT("Vectors.size()!=Duals.size()");
    size_t testSize=std::min(size_t(sqrt(1.e8/Vectors.size())),Vectors.size());
    std::vector<size_t> iSub;
    for(size_t k=0;k<Vectors.size();k++)iSub.push_back(k);
    while(iSub.size()>testSize)
        iSub.erase(iSub.begin()+rand()%iSub.size());
    for(auto i: iSub)
        for(auto j: iSub)
            if(std::abs(Duals[i]->innerProduct(Vectors[j],true)-double(i==j))>Epsilon)
                DEVABORT(Sstr+"Duals do not match Vectors by"+std::abs(Duals[i]->innerProduct(Vectors[j],true)-double(i==j)));

    if(iSub.size()<Vectors.size())PrintOutput::DEVmessage(Sstr+"Tested Duals, random  sample of"+iSub.size()+"out of total"+Vectors.size());
    else                          PrintOutput::DEVmessage(Sstr+"Tested Duals"+Vectors.size());
    PrintOutput::outputLevel("restore");
}

void ProjectSubspace::_construct(std::vector<const Coefficients *> Vectors, std::vector<const Coefficients *> Duals,size_t BeginOrthonormal)
{
    if(Vectors.size()!=Duals.size())DEVABORT(Sstr+"unequal number of Vectors and Duals"+Vectors.size()+Duals.size());
    name="Project"+tools::str(Vectors.size());
    if(Vectors.size()==0)return;
    if(BeginOrthonormal*Vectors.size()>1.e6)PrintOutput::DEVmessage(Sstr+"setting up large projector from"+BeginOrthonormal+"non-ON vectors");
    iIndex=Vectors[0]->idx();
    jIndex=iIndex;

    // check duals
    if(Vectors[0]->idx()!=Duals[0]->idx())DEVABORT(Sstr+"unequal Index on Vectors and Duals"+Vectors[0]->idx()+Duals[0]->idx());
    checkDuals({Vectors.begin()+BeginOrthonormal,Vectors.end()},{Duals.begin()+BeginOrthonormal,Duals.end()},1.e-10);

    std::vector<std::vector<std::vector<const Coefficients*> > > gVecs=groupByNonzeros(Vectors,Duals,_sorting);

    // build a grouped subspace index
    Index* idx=new Index(std::vector<const BasisAbstract*>(1,BasisAbstract::factory("Vector:"+tools::str(gVecs[0].size()))),{"Subspace"});
    for(size_t k=0;k<idx->basis()->size();k++)
        idx->childReplace(k,new Index(std::vector<const BasisAbstract*>(1,BasisAbstract::factory("Vector:"+tools::str(gVecs[0][k].size()))),{"SubSubspace"}));
    std::vector<std::string> dum;
    idx->setFloorAuto(dum);
    idx->sizeCompute();

    // for now, we consider projections as global (can be improved)
    _subspaceIndex=idx;
    _subspaceC.reset(new Coefficients(_subspaceIndex));

    _mapTo.reset  (new OperatorTree("mapToFullFromContracted",iIndex,_subspaceIndex));
    _mapFrom.reset(new OperatorTree("mapFromFullToContracted",_subspaceIndex,iIndex));

    for(size_t k=0;k<gVecs[0].size();k++){
        _mapTo->childAdd(new OperatorTree("block",iIndex,_subspaceIndex->child(k)));
        _mapFrom->childAdd(new OperatorTree("block",_subspaceIndex->child(k),iIndex));
        buildMap(_mapTo->childBack(),gVecs[0][k]);
        buildMap(_mapFrom->childBack(),gVecs[1][k]);
    }
    _mapTo->purge(1.e-12);
    _mapFrom->purge(1.e-12);

    if(BeginOrthonormal>0){
        std::vector<Eigen::Triplet<std::complex<double>>> list;
        for(size_t j=0;j<BeginOrthonormal;j++){
            for(size_t i=0;i<Duals.size();i++){
                list.push_back(Eigen::Triplet<std::complex<double> >(i,j,Duals[i]->innerProduct(Vectors[j],true)));
                if(i>=BeginOrthonormal)
                    list.push_back(Eigen::Triplet<std::complex<double> >(j,i,Duals[j]->innerProduct(Vectors[i],true)));
            }
        }
        // fill up diagonal
        for(size_t k=BeginOrthonormal;k<Vectors.size();k++)
            list.push_back(Eigen::Triplet<std::complex<double> >(k,k,1));

        Eigen::SparseMatrix<std::complex<double>> sov(Duals.size(),Vectors.size());
        sov.setFromTriplets(list.begin(),list.end());

        // permute the sparse matrix into _sorting
        Eigen::PermutationMatrix<Eigen::Dynamic> perm(sov.rows());
        for(size_t k=0;k<_sorting.size();k++)perm.indices()[_sorting[k]]=k;
        sov=sov.twistedBy(perm);

        // get the sparse LU decomposition
        _lu.reset(new Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>,Eigen::COLAMDOrdering<int>>(sov));

        // need vector collected on master

        for(Index * jdx=const_cast<Index*>(_subspaceIndex->firstFloor());jdx;jdx=jdx->nodeNext())
            Parallel::setIndexOwner(jdx,MPIwrapper::master());
    }

    ParallelOperator::setDistribution(_mapTo.get());
    ParallelOperator::setDistribution(_mapFrom.get());

    if(not verify())DEVABORT("basic projector does not work");
}

void ProjectSubspace::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    if(not _subspaceC)const_cast<ProjectSubspace*>(this)->_subspaceC.reset(new Coefficients(_subspaceIndex));
    _mapFrom->apply(A,Vec,0.,*_subspaceC);
    if(_lu!=0){
        if(MPIwrapper::Size()>1)DEVABORT("ProjectSubspace::apply is not parallel yet");
        if((int)_subspaceC->size()!=_lu->rows())DEVABORT("sizes do not match");
        Eigen::Map<Eigen::VectorXcd>(_subspaceC->data(),_subspaceC->size())
                =_lu->solve(Eigen::Map<Eigen::VectorXcd>(_subspaceC->data(),_subspaceC->size()));
    }
    _mapTo->apply(1.,*_subspaceC,B,Y);
}

void ProjectSubspace::applyDual(std::complex<double> A, const Coefficients &Dual, std::complex<double> B, Coefficients &Y) const{
    // rarely used - if needed, implement more efficiently
    if(not _subspaceC)const_cast<ProjectSubspace*>(this)->_subspaceC.reset(new Coefficients(_subspaceIndex));
    Dual.idx()->inverseOverlap()->apply(1.,Dual,0,*Dual.idx()->inverseOverlap()->tempLHS());
    _mapFrom->apply(A,*Dual.idx()->inverseOverlap()->tempLHS(),0.,*_subspaceC.get());
    if(_lu!=0){
        if(MPIwrapper::Size()>0)DEVABORT("ProjectSubspace::apply is not parallel yet");
        if((int)_subspaceC->size()!=_lu->rows())DEVABORT("sizes do not match");
        Eigen::Map<Eigen::VectorXcd>(_subspaceC->data(),_subspaceC->size())
                =_lu->solve(Eigen::Map<Eigen::VectorXcd>(_subspaceC->data(),_subspaceC->size()));
    }
    _mapTo->apply(1.,*_subspaceC.get(),0,*Dual.idx()->inverseOverlap()->tempLHS());
    Dual.idx()->overlap()->apply(1.,*Dual.idx()->inverseOverlap()->tempLHS(),B,Y);
}

std::vector<Coefficients> ProjectSubspace::orbitals(std::vector<int> Select){
    std::vector<int>select(Select);
    if(select.size()==0)
        for(int k=0;k<dim();k++)select.push_back(k);

    std::vector<Coefficients>result;
    for(int k: select){
        _subspaceC->setToZero();
        _subspaceC->orderedData()[k]=1.;
        result.push_back(Coefficients(iIndex));
        _mapTo->apply(1.,*_subspaceC,0.,result.back());
    }
    return result;
}

bool ProjectSubspace::verify() const {
    if(MPIwrapper::Size()>1){
        PrintOutput::DEVwarning("ProjectSubspace::verify() is not parallel - skipped");
        return true;
    }
    Coefficients c(iIndex),d(iIndex);
    apply(1.,c,0.,d);
    apply(1.,d,0.,c);
    if(not (d-=c).isZero(1.e-12))
        PrintOutput::DEVwarning(Sstr+"spectral projectors failed for "+name+"by"+d.cMaxNorm());
    return d.isZero(1.e-9);
}


ProjectSubspace::~ProjectSubspace(){
}
