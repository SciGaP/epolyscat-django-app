// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "eigenSolver.h"
#include "eigenInverseIter.h"

#include "operatorTree.h"
#include "index.h"
#include "coefficients.h"
#include "coefficientsSparse.h"
#include "str.h"
#include "printOutput.h"
#include "inverse.h"
#include "operatorAbstract.h"
#include "parallel.h"
#include "parallelOperator.h"
#include "mpiWrapper.h"

#include "eigenTools.h"
#include "timer.h"
#include "timeCritical.h"

#include "resolvent.h"
#include "algebra.h"

//debug
#include "coefficients.h"
#include "overlapDVR.h"

#include "indexNew.h"
#include "matrixBlocking.h"
#include "matrixEigenBlock.h"

TIMER(eigen0,)
TIMER(eigen1,)
TIMER(eigen2,)
TIMER(eigen3,)
TIMER(eigen4,)

using namespace std;

EigenSolver::EigenSolver():EigenSolverAbstract(),_method("auto"),_fullVectors(false){withSelect("All");}

EigenSolver::EigenSolver(double Emin, double Emax, int Nmax, bool RightVectors, bool DualVectors, bool ExcludeRange, string Method)
    :EigenSolver()
{
    _method=Method;
    makeSelect(Nmax,Emin,Emax,ExcludeRange);
    computeLeftVectors(DualVectors);
    computeRightVectors(RightVectors);
}

EigenSolver::~EigenSolver(){};

void EigenSolver::_compute()
{
    LOG_PUSH("compB1");
    std::string log=_op->name+tools::str(depth());

    // descend in block-diagonal part
    const OperatorTree* OpTree=dynamic_cast<const OperatorTree*>(_op);
    const OperatorTree* OvrTree=dynamic_cast<const OperatorTree*>(_ovr);
    // re-distribute by diagonal blocks
    ParallelOperator parH(OpTree),parO(OvrTree);
    ParallelOperator::SavedHost savedH(&parH),savedO(&parO);
    parH.reDistribute(ParallelOperator::DiagonalBlockHost(OpTree));
    parO.reDistribute(ParallelOperator::DiagonalBlockHost(OvrTree));
    LOG_POP();
    LOG_PUSH("compB2");

    // diagonalize block-wise
    _computeBlock(OpTree,OvrTree);
    LOG_POP();
    LOG_PUSH("compB3");

    // restore original distribution
    parH.reDistribute(savedH);
    parO.reDistribute(savedO);

    //CAUTION: this may become very (unnecessarily) large
    if(MPIwrapper::Size()>1 and ParallelOperator::getHost(_op)==ParallelOperator::distributed){
        for(EigenSolver* leaf=firstLeaf();leaf!=0;leaf=leaf->nextLeaf()){
            Parallel::gatherAllEigen(leaf->_op->iIndex,leaf->_eigenvalues,leaf->_rightVectors,leaf->_dualVectors);
        }
    }
    LOG_POP();
    LOG_PUSH("compB4");
    select(_select);
    LOG_POP();
}

// find overlap block from index hierarchy
const OperatorAbstract* ovrBlock(const Index* iIndex,const Index* jIndex){
    const Index* iRoot=iIndex;
    for(;iRoot and iRoot->overlap()==0;iRoot=iRoot->parent());
    if(iRoot){
        for(const OperatorTree* o=dynamic_cast<const OperatorTree*>(iRoot->overlap());o!=0;o=o->nodeNext()){
            if(o->iIndex==iIndex and o->jIndex==jIndex)return o;
        }
    }
    DEVABORT("no overlap block found");
    return 0;
}



void EigenSolver::_computeBlock(const OperatorTree* OpTree, const OperatorAbstract *Ov)
{
    _ovr=Ov;
    _op=OpTree;
    if(!_ovr){
        IndexOverlap::set(_op->iIndex);
        _ovr=IndexOverlap::get(_op->iIndex);
    }
    if(parent() and _op->iIndex->size()<1000)PrintOutput::outputLevel("restore");
    const OperatorTree* ovTree=dynamic_cast<const OperatorTree*>(Ov);
    if(ovTree!=0 and OpTree->isBlockDiagonal() and ovTree->isBlockDiagonal()){
        for(size_t k=0;k<OpTree->childSize();k++){
            for(size_t l=0;l<k;l++){
                if(OpTree->child(k)->iIndex==OpTree->child(l)->iIndex)
                    DEVABORT("not implemented for equal index blocks\n"+OpTree->str());
            }

            childAdd(new EigenSolver(eMin(),eMax(),numberEigenv(),_computeRightVectors,_computeLeftVectors,excludeRange(),_method));
            childBack()->withSelect(_select);
            childBack()->_computeBlock(OpTree->child(k),OpTree->child(k)->iIndex->overlap());
        }
    }
    else {

        if(parent()==0)
            PrintOutput::DEVmessage(Str("Eigensolver for problem size"," ")+OpTree->iIndex->sizeStored()+"("+SEP("")+autoMethod(_method)+")");
        else if(nSibling()==0){
            vector<int>siz;
            for(const OperatorTree*o=OpTree;o!=0;o=o->nodeRight())siz.push_back(o->iIndex->sizeStored());
        }

        int host = ParallelOperator::getHost(OpTree);

        if(host==MPIwrapper::Rank() or host==ParallelOperator::all)
        {
            if(OpTree->iIndex->sizeStored()>1 and autoMethod(_method).find("Arpack")==0){


                if(OpTree->iIndex->size()>200 and numberEigenv()*2>int(OpTree->iIndex->size()))
                    ABORT(Sstr+"Too many eigenvaluse for Arpack solve \n"
                               "requested: "+numberEigenv()+" problem size:"+OpTree->iIndex->size());

                // number of vectors
                int nVec=std::min(numberEigenv(),int(OpTree->iIndex->size())-1);

                PrintOutput::message(Sstr+"using"+autoMethod(_method)+"for"+nVec+"eigenvalues out of"+OpTree->iIndex->size());

                string which="SmallReal";
                if(eMin()>0. and eMax()>=DBL_MAX/4.)which="LargeAbs";
                vector<complex<double> > eval;
                vector<Coefficients*> evec;
                bool crit=timeCritical::coefficients;
                timeCritical::coefficients=false;

                if(autoMethod(_method)=="Arpack"){
                    Arp a(OpTree);
                    a.eigen(eval,evec,nVec,which,false);
                }
                else if(autoMethod(_method)=="ArpackShiftInverse"){
                    std::complex<double> eguess=0.5*(eMin()+eMax());
                    if(_select.find("Nearest[")==0){
                        std::vector<std::string> sel=tools::splitString(tools::stringInBetween(_select,"[","]"),',');
                        eguess=std::complex<double>(Algebra(sel[1]).val(0).real(),Algebra(sel[2]).val(0).real());
                    }
                    ArpShiftInverse a(OpTree,eguess);
                    a.eigen(eval,evec,nVec,which,false);
                }
                timeCritical::coefficients=crit;

                // select desired energy range
                for(size_t k=0;k<eval.size();k++){
                    if(eMin()<=eval[k].real() and eval[k].real()<=eMax()){
                        _eigenvalues.push_back(eval[k]);
                        if(evec.size()>k)_rightVectors.push_back(evec[k]);
                    }
                    else
                        if(evec.size()>k)delete evec[k];
                }

                // get dual basis
                symmetricDuals();

            }

            else if (OpTree->iIndex->sizeStored()<=1 or autoMethod(_method)=="Lapack"){
                if(!_ovr)DEVABORT("empty overlap for "+OpTree->iIndex->hierarchy());

                LOG_PUSH("eigenA1");
                Eigen::MatrixXcd eVal;
                Eigen::SparseMatrix<std::complex<double>> rVec,dVec;
                Eigen::SparseMatrix<std::complex<double>> sOp=OpTree->matrixSparse(true);
                Eigen::SparseMatrix<std::complex<double>> sOv=_ovr->matrixSparse(true);
                LOG_POP();

                if(not (_computeLeftVectors or _computeRightVectors)){
                    tRecX::eigenBlock<Eigen::SparseMatrix<std::complex<double>>>(sOp,eVal,false,rVec,false,dVec,sOv);
                    for(int k=0;k<eVal.size();k++)
                        if(eMin()<=eVal.data()[k].real() and eVal.data()[k].real()<=eMax())
                            _eigenvalues.push_back(eVal.data()[k]);
                }
                else {
                    LOG_PUSH("eigenA2");
                    // do not assume simple relation between left and right eigenvectors
                    tRecX::eigenBlock<Eigen::SparseMatrix<std::complex<double>>>(sOp,eVal,true,rVec,true,dVec,sOv);
                    LOG_POP();
                    LOG_PUSH("eigenA3");
                    expandInto(_rightVectors,_eigenvalues,eVal.data(),rVec,eMin(),eMax(),excludeRange());
                    expandInto( _dualVectors,_eigenvalues,eVal.data(),dVec,eMin(),eMax(),excludeRange());
                    LOG_POP();
                }
            }

            else if (autoMethod(_method)=="InverseIteration"){
                EigenInverseIter invit;
                invit.withGuesseigenvalue(guessEigenvalue());
                PrintOutput::message("using "+autoMethod(_method)+" with selection "+invit.selection());
                invit.compute(_op,_ovr);
                _eigenvalues=invit.eigenvalues();
                _rightVectors.assign(1,new Coefficients(*invit.rightVectors()[0]));
                _dualVectors.assign(1,new Coefficients(_op->iIndex));
                _ovr->apply(1.,*_rightVectors[0],0.,*_dualVectors[0]);
            }
            else ABORT("unknown EigenSolver method: "+_method);
            select(_select);

            LOG_PUSH("eigenA4");
            if(_computeRightVectors){
                orthonormalizeDegenerate(1.e-2,_eigenvalues,_dualVectors,_rightVectors);
                if(not verify())PrintOutput::DEVwarning("eigensolving failed");
            }
            LOG_POP();
        }
        select(_select);
    }
    PrintOutput::outputLevel("low");
}

void EigenSolver::collect(){
    if(_eigenvalues.size())return; // already collected
    int dep=_op->iIndex->depth(); // determines level on which to collect (usually top)
    _epsVerified=0.;
    for(EigenSolver * leaf=firstLeaf();leaf!=0;leaf=leaf->nextLeaf(this)){
        _epsVerified=std::max(_epsVerified,leaf->_epsVerified);
        Coefficients cf(_op->iIndex);
        for(size_t k=0;k<leaf->_eigenvalues.size();k++){
            _eigenvalues.push_back(leaf->_eigenvalues[k]);

            if(_computeRightVectors){
                vector<unsigned int>lidx(leaf->_op->iIndex->index()); // multi-index of compute block
                lidx.erase(lidx.begin(),lidx.begin()+dep);            // remove top part of multi-index
                if(_fullVectors){
                    _rightVectors.push_back(new Coefficients(_op->iIndex));
                    *_rightVectors.back()->nodeAt(lidx)=*leaf->_rightVectors[k];

                    _dualVectors.push_back(new Coefficients(_op->iIndex));
                    *_dualVectors.back()->nodeAt(lidx)=*leaf->_dualVectors[k];
                } else {
                    Coefficients* cNode=cf.nodeAt(lidx);
                    *cNode=*leaf->_rightVectors[k];
                    _rightVectors.push_back(new CoefficientsSparse(cf));
                    *cNode=*leaf->_dualVectors[k];
                    _dualVectors.push_back(new CoefficientsSparse(cf));
                }
            }
        }
    }
    if(_epsVerified==0.)_epsVerified=DBL_MAX;
    phaseFix();
}
EigenSolver & EigenSolver::fullVectors(){
    PrintOutput::DEVwarning("constructing full eigenvectors - hack only for DiscretizationSpectral - fix at one point");
    _fullVectors=true;
    return *this;
}

std::vector<Coefficients*> EigenSolver::rightVectors() {
    if(!_rightVectors.size())collect();
    return _rightVectors;
}
std::vector<Coefficients*> EigenSolver::dualVectors() {
    if(!_dualVectors.size())collect();
    return  _dualVectors;
}

vector<complex<double> > EigenSolver::eigenvalues() {
    if(!_eigenvalues.size())collect();
    return _eigenvalues;
}


// Eigen::SparseMatrix version of Index::unGlobal
static void unGlobal(const Index* Idx, const Eigen::SparseMatrix<std::complex<double>> &Eigen, vector<Coefficients*> & Evec, vector<int> cols) {
    vector<unsigned int>cI(Idx->contractedNumbering());
    for (unsigned int col=0; col<cols.size(); col++) {
        Evec.push_back(new Coefficients(Idx));
        const Coefficients * leaf=Evec.back()->firstLeaf();
        while(leaf!=0){
            unsigned int pos=leaf->idx()->posIndex(Idx);
            for(unsigned int k=0;k<leaf->size();k++){
                leaf->floorData()[k]=Eigen.coeff(cI[pos+k],cols[col]);
            }
            leaf=leaf->nextLeaf();
        }
        // contracted dual vectors carry the contraced overlap matrix
        // remove that extra factor 2 from margins
        Evec.back()->makeContinuous(sqrt(0.5));
    }
}

void EigenSolver::expandInto(std::vector<Coefficients *> &CVec, std::vector<std::complex<double> > &Val,
                             const std::complex<double> *PVal,  const Eigen::SparseMatrix<std::complex<double>> & Vec, double Emin, double Emax,
                             bool ExcludeRange)
{
    // list of selected eigenvalues
    Val.clear();
    vector<int> cols;
    for(int k=0;k<Vec.cols();k++)
        if((Emin<=PVal[k].real() and PVal[k].real()<=Emax) != ExcludeRange){
            cols.push_back(k);
            Val.push_back(PVal[k]);
        }

    // expand selected values
    unGlobal(_op->iIndex,Vec,CVec,cols);

}

complex<double> EigenSolver::Arp::shift=0;

EigenSolver::Arp::Arp(const OperatorAbstract *O)
    :Arpack(1.e-9,5000,true),o(O),x(Coefficients(O->jIndex)),y(Coefficients(O->iIndex)),tmp(Coefficients(O->iIndex)){
    x.pointerToC(pX);
    y.pointerToC(pY);
    tmp.setToZero(); // make sure =0
    y.setToZero(); // make sure Y=0
    _lvec=y.size();
}

void EigenSolver::Arp::eigen(vector<std::complex<double> > &Eval, vector<Coefficients *> &Rvec,
                             unsigned int Nvec, const string &Which, bool Restart){

    complex<double> defaultShift=EigenSolver::Arp::shift;
    if(Restart){
        // initialize with starting vector
        if(Rvec.size()!=1)ABORT("for restarting, supply exactly one Rvec as starting vector");
        y=*Rvec[0];
        //delete Rvec[0];
        Rvec.clear();
        rvec.resize(_lvec);
        for(unsigned int i=0;i<_lvec;i++)rvec[i]=*(pY[i]);
    }
    else if(Rvec.size()>0)
        ABORT("must enter eigen with empty Rvec, unless Restart");

    if(Which=="SmallAbs"){
        if(Eval.size()!=1)ABORT("for Which=SmallAbs provide guess energy value in Eval");
        EigenSolver::Arp::shift=Eval[0];
    }
    else if(Which=="SmallReal"){
        // in FEM, there may be 0-values due to projection
        EigenSolver::Arp::shift=100.; // if true ground statet is >0, need to shift safely(?) to negative values
    }

    // let ArpackFunction compute the vectors
    Arpack::eigenIter(Nvec,Which,Restart,true);

    // move into Rvec
    Eval.resize(Nvec);
    for(unsigned int k=0;k<Nvec;k++){
        Eval[k]=eval[k]+EigenSolver::Arp::shift;

        for(unsigned int i=0;i<_lvec;i++)*(pY[i])=rvec[i+k*_lvec];
        Rvec.push_back(new Coefficients(y));
    }

    // reset in case is had been changed
    EigenSolver::Arp::shift=defaultShift;
}


void EigenSolver::Arp::apply(const std::complex<double> *X, std::complex<double> *Y) {

    // copy to Coefficients
    for (unsigned int i=0;i<pX.size();i++)*(pX[i])=*(X+i);

    // y=S^-1 Op x - shift*x
    x.makeContinuous();
    o->apply(1.,x,0.,tmp);
    tmp.makeContinuous();
    if(not o->iIndex->inverseOverlap()){
        PrintOutput::DEVwarning("lost overlap/inverse overlap - this should be fixed");
        dynamic_cast<IndexNew*>(const_cast<Index*>(o->iIndex))->buildOverlap();
    }

    o->iIndex->inverseOverlap()->apply(1.,tmp,0.,y);
    y.axpy(-EigenSolver::Arp::shift,&x);

    // copy back and reset y to zero
    for (unsigned int i=0;i<pX.size();i++)*(Y+i)=*(pY[i]);
}

EigenSolver::ArpShiftInverse::ArpShiftInverse(const OperatorAbstract* Op,const std::complex<double> Eguess):Arp(Op){
    const OperatorTree* op=dynamic_cast<const OperatorTree*>(Op);
    if(op)_resolv.reset(new Resolvent(op,Eguess));
    else DEVABORT("ArpShifInverse only for OperatorTree");
}

void EigenSolver::ArpShiftInverse::apply(const std::complex<double> *X, std::complex<double> *Y) {

    // copy to Coefficients
    for (unsigned int i=0;i<pX.size();i++)*(pX[i])=*(X+i);

    // y=(Op-Eguess)^-1 S x
    x.makeContinuous();
    o->iIndex->overlap()->apply(1.,x,0.,tmp);
    tmp.makeContinuous();
    _resolv->apply(1.,tmp,0.,y);

    // copy back
    for (unsigned int i=0;i<pX.size();i++)*(Y+i)=*(pY[i]);
}

void EigenSolver::ArpShiftInverse::eigen(std::vector<std::complex<double> > &Eval, std::vector<Coefficients *> &Rvec, unsigned int Nvec, const std::string &Which, bool Restart){

    if(Restart){
        // initialize with starting vector
        if(Rvec.size()!=1)ABORT("for restarting, supply exactly one Rvec as starting vector");
        y=*Rvec[0];
        //delete Rvec[0];
        Rvec.clear();
        rvec.resize(_lvec);
        for(unsigned int i=0;i<_lvec;i++)rvec[i]=*(pY[i]);
    }
    else if(Rvec.size()>0)
        ABORT("must enter eigen with empty Rvec, unless Restart");

    // let ArpackFunction compute the vectors
    Arpack::eigenIter(Nvec,Which,Restart,true);

    // back-transform spectrum and move into Rvec
    Eval.resize(Nvec);
    for(unsigned int k=0;k<Nvec;k++){
        // exclude zeros that may arrise from continuity conditions
        if(abs(eval[k])>1.e-12){
            Eval[k]=_resolv->z()+1./eval[k];
            for(unsigned int i=0;i<_lvec;i++)*(pY[i])=rvec[i+k*_lvec];
            Rvec.push_back(new Coefficients(y));
        }
    }
}

