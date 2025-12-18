// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "derivativeFlat.h"
#include "mpiWrapper.h"

#include "printOutput.h"
#include "str.h"
#include "parallel.h"

#include "operatorFloor.h"
#include "operatorSingle.h"
#include "index.h"
#include "coefficientsFloor.h"
#include "discretizationSpectral.h"
#include "discretizationSpectralProduct.h"
#include "operatorFloor.h"
#include "operatorAbstractProduct.h"
#include "parameters.h"
#include "operatorDiagonal.h"
#include "readInput.h"
#include "derivativeBlock.h"
#include "parallelProcess.h"
#include "coefficients.h"
#include "coefficientsGlobal.h"
#include "parallelOperator.h"
#include "projectSubspace.h"

#include "log.h"
#include "operatorMeanEE.h"
#include "operatorNonLin.h"

#include "mpiWrapper.h"
#include "asciiFile.h"

#include "inverseDvr.h"
#include "toolsPrint.h"
#include "timeCritical.h"

//#include "operatorMeanEE.h"


using namespace std;
bool DerivativeFlat::applyFlat=false;
complex<double> DerivativeFlat::Arp::shift=complex<double>(0.);
complex<double> DerivativeFlat::Arp::shiftZero=complex<double>(1000.);

DerivativeFlat::~DerivativeFlat(){
    delete op;
    delete setupTemp;
    delete setupXY;
}
DerivativeFlat::FlattenedOperatorTree::~FlattenedOperatorTree(){delete par;}

const Coefficients & DerivativeFlat::lhsVector() const {return *setupXY;}

// make a list of operators to update
//std::set<OperatorTree*> _nonLinOps;
//void DerivativeFlat::setNonlinUpdate(){
//    if(op->par->_process.size()>1)ABORT("OperatorNonLin cannot be used in parallel");
//    OperatorFloorNonLin* floorNL;
//    for(ParallelGrain* g: op->par->grain)
//        for(DerivativeBlock* b: g->block){
//            OperatorTree* nonLin=const_cast<OperatorTree*>(b->oLeaf);
//            if(floorNL=dynamic_cast<OperatorFloorNonLin*>(nonLin->floor())){
//                while(nonLin and not dynamic_cast<OperatorNonLin*>(nonLin))nonLin=nonLin->parent();
//                if(nonLin)_nonLinOps.insert(nonLin);
//                else _nonLinOps.insert(const_cast<OperatorTree*>(b->oLeaf));
//        }
//    }
//}


static bool updateFirst=true;
static bool updateOff=false;
void DerivativeFlat::updateNonLin(double Time,Coefficients *C){

    if(updateFirst){
        DEVABORT("check logics");
        if(updateOff)PrintOutput::warning("mean field update off");
        updateFirst=false;
    }
    if(updateOff)return;

    OperatorFloorNonLin* floorNL;
    for(ParallelGrain* g: op->par->grain)
        for(DerivativeBlock* b: g->block){
            if((floorNL=dynamic_cast<OperatorFloorNonLin*>(const_cast<OperatorTree*>(b->oLeaf)->floor()))){
                if(op->par->_process.size()>1)ABORT("OperatorNonLin cannot be used in parallel");
                floorNL->updateNonLin(Time,C);
            }
        }
}


void DerivativeFlat::update(double Time, const Coefficients* CurrentVec)
{
    _time=Time;
    Parameters::update(Time,false);
    //if(CurrentVec)updateNonLin(Time,const_cast<Coefficients*>(CurrentVec));
}

void DerivativeFlat::test(){
    // note: the test has only be set up w/o projection
    if(o->isHuge())return;
    if(projection!=0)return;

    Parameters::updateToOne();
    UseMatrix dMat,oMat,sMat;
    o->matrix(oMat);
    o->iIndex->inverseOverlap()->matrix(sMat);
    oMat=sMat*oMat;
    oMat*=complex<double>(0,-1.);
    matrix(dMat);
    if(MPIwrapper::isMaster()){
        if(not (dMat-oMat).isZero(1.e-10)){
            dMat.print("dmat",2);
            oMat.print("omat",2);
            (oMat-=dMat).print("diff",0);
            ABORT(Str("test failed with error")+oMat.maxAbsVal());
        }
        else cout<<"OK "<<endl;
    }
    Parameters::restoreToTime();
}


void DerivativeFlat::FlattenedOperatorTree::add(complex<double>* Fac, const OperatorTree* Op,
                                                Coefficients *ICoeff, Coefficients *JCoeff, bool AbsorbInverse, double ApplyEpsilon, const vector<unsigned int> &ISort,
                                                const vector<unsigned int> &JSort){

    if(Op->floor()!=0){
        vector<unsigned int> blockSort;
        for(unsigned int k=0;k<JSort.size();k++)blockSort.push_back(JCoeff->idx()->index()[JSort[k]]);
        for(unsigned int k=0;k<ISort.size();k++)blockSort.push_back(ICoeff->idx()->index()[ISort[k]]);

        double eps=DBL_MAX;
        if(Op->floor()->norm()!=0)eps=ApplyEpsilon/Op->floor()->norm();
        double nrm=JCoeff->norm(); // pointer does not survive - must not be used!
        blocks.push_back(DerivativeBlock(blockSort, Op, ICoeff, JCoeff, eps, &nrm));
        if(blocks.back().eps==DBL_MAX){
            PrintOutput::DEVwarning(Str("zero-block in"," ")+Op->root()->name+"- should have been purged\n"
                                    +Op->name
                                    +Op->iIndex->index()+"|"+Op->jIndex->index()
                                    +Op->iIndex->basis()->str()
                                    +Op->jIndex->basis()->str()
                                    +"["+MPIwrapper::Rank()+"]");
        }
    }
    if(Op->isLeaf() and Op->floor()==0){
        MPout<<(Str("empty operator ")+Op->str()+Op->iIndex->index()+Op->jIndex->index()+"\n");
    }

    // descend in hierarchy
    for(unsigned int k=0;k<Op->childSize();k++){
        Coefficients*iC=ICoeff,*jC=JCoeff;

        if(Op->iIndex!=Op->child(k)->iIndex and Op->iIndex->childSize()>0)
            iC=ICoeff->child(Op->iIndex->nSub(Op->child(k)->iIndex));
        if(Op->jIndex!=Op->child(k)->jIndex and Op->jIndex->childSize()>0)
            jC=JCoeff->child(Op->jIndex->nSub(Op->child(k)->jIndex));
        add(Fac, Op->child(k), iC, jC, AbsorbInverse, ApplyEpsilon, ISort, JSort);
    }

    // after all blocks added, sort
    // using "stable_sort" as "sort" did not function properly in some cases
    // cause for this behavior is unclarified
    if(Op->iIndex->isRoot() and Op->jIndex->isRoot())
        std::stable_sort(blocks.begin(), blocks.end(), DerivativeBlock::lessEqual);
}

DerivativeFlat::FlattenedOperatorTree::FlattenedOperatorTree(const OperatorAbstract* Map, bool AbsorbInverse,
                                                             double ApplyEpsilon, Coefficients* ICoeff, Coefficients* JCoeff, string SendRecv){

    if(Map->idx()!=Map->jdx()){
        if(SendRecv!="send" and SendRecv!="receive" and SendRecv!="either")
            DEVABORT("idx()!=jdx(), need to use send or receive, have: "+SendRecv);
    }
    const OperatorTree* treeMap = dynamic_cast<const OperatorTree*>(Map);

    // make sure we are sync'd
//    ParallelOperator::sync(const_cast<OperatorTree*>(treeMap));

    if(treeMap==0)ABORT("need spectral map as operators tree");
    add(0, treeMap, ICoeff, JCoeff, AbsorbInverse, ApplyEpsilon,
        Parallel::indexSort[ICoeff->idx()], Parallel::indexSort[JCoeff->idx()]);

    par=new Parallel(MPIwrapper::Size(),3);
    par->addBlocks(blocks,SendRecv);
    if(MPIwrapper::Size()>1)PrintOutput::verbatim("\nParallel layout of "+Map->name+"\n"+par->str()+"\n");
    if(par->str().find("(unused)")!=string::npos)PrintOutput::warning("unused process in"+treeMap->name);

    ParallelOperator parOp(treeMap);
    parOp.reDistribute(ParallelOperator::ProcessHost(par->_process));
}

DerivativeFlat::ProjectionSingle::ProjectionSingle(const DiscretizationSpectral* ProjectionDisc, Coefficients* Coeff,
                                                   double ApplyEpsilon)
    :cSpec(new Coefficients(ProjectionDisc->mapFromParent()->idx()))
    ,cXY(Coeff),cTemp(new Coefficients(Coeff->idx()))
{

    globalXY = CoefficientsGlobal::view(cXY);
    localXY = CoefficientsLocal::view(globalXY);

    globalTemp = CoefficientsGlobal::view(cTemp);
    localTemp = CoefficientsLocal::view(globalTemp);

    mapFrom = new FlattenedOperatorTree(ProjectionDisc->mapFromParent(), false, ApplyEpsilon, cSpec, cXY, "send");
    mapTo = new FlattenedOperatorTree(ProjectionDisc->mapToParent(), false, ApplyEpsilon, cXY, cSpec, "receive");

    globalSpec=CoefficientsGlobal::view(cSpec);
    localSpec=CoefficientsLocal::view(globalSpec);

}

DerivativeFlat::ProjectionSingle::ProjectionSingle(const ProjectSubspace* Projection, Coefficients* Coeff,
                                                   double ApplyEpsilon)
    :cSpec(new Coefficients(Projection->mapFrom()->idx())),cXY(Coeff),cTemp(new Coefficients(Coeff->idx())),_lu(0)
{
    cSpec->setToZero();

    globalXY = CoefficientsGlobal::view(cXY);
    localXY = CoefficientsLocal::view(globalXY);

    globalTemp = CoefficientsGlobal::view(cTemp);
    localTemp = CoefficientsLocal::view(globalTemp);

    mapFrom = new FlattenedOperatorTree(Projection->mapFrom().get(), false, ApplyEpsilon, cSpec, cXY, "either");
    mapTo =   new FlattenedOperatorTree(Projection->mapTo().get(),   false, ApplyEpsilon, cTemp,cSpec, "either");
    _lu=Projection->lu();

    globalSpec=CoefficientsGlobal::view(cSpec);
    localSpec=CoefficientsLocal::view(globalSpec);
}

TIMER(proApplyA,)
TIMER(proApplyB,)
TIMER(proApplyC,)
void DerivativeFlat::ProjectionSingle::apply(){
    if(cSpec->size()==0)return;
    cSpec->setToZero();
    localSpec->setToZero();
    localTemp->setToZero();
    globalSpec->setToZero();
    globalTemp->setToZero();

    mapFrom->par->apply(1.);
    if(_lu!=0){
        if(_lu->cols()!= (int) cSpec->size())DEVABORT(Sstr+"BAD"+_lu->rows()+_lu->cols()+cSpec->size()
                                                          +cSpec->idx()->str(Tree_defaultKind,3));
        Eigen::Map<Eigen::VectorXcd>(cSpec->data(),cSpec->size())
                =_lu->solve(Eigen::Map<Eigen::VectorXcd>(cSpec->data(),cSpec->size()));
    }
    mapTo->par->apply(1.);

    *localXY -= *localTemp;
}


DerivativeFlat::ProjectionProduct::ProjectionProduct(const DiscretizationSpectralProduct* ProjectionDisc, Coefficients* Coeff,
                                                     double ApplyEpsilon):
    cXY(Coeff),
    cTemp(new Coefficients(Coeff->idx())){

    globalXY = CoefficientsGlobal::view(cXY);
    localXY = CoefficientsLocal::view(globalXY);

    globalTemp = CoefficientsGlobal::view(cTemp);
    localTemp = CoefficientsLocal::view(globalTemp);

    for(size_t i=0; i<ProjectionDisc->factors.size(); i++){
        // determine index owner for spectral blocks (primitive for now)

        cSpec.push_back(new Coefficients(ProjectionDisc->factors[i]->mapFromParent()->iIndex));
        mapFrom.push_back(
                    new FlattenedOperatorTree(
                        ProjectionDisc->factors[i]->mapFromParent(),
                        false,
                        ApplyEpsilon,
                        cSpec[i],
                        cXY,
                        "either"
                        )
                    );
        mapTo.push_back(
                    new FlattenedOperatorTree(
                        ProjectionDisc->factors[i]->mapToParent(),
                        false,
                        ApplyEpsilon,
                        cTemp,
                        cSpec[i],
                        "either"
                        )
                    );

        globalSpec.push_back(CoefficientsGlobal::view(cSpec.back()));
        localSpec.push_back(CoefficientsLocal::view(globalSpec.back()));
    }
}


TIMER(applyProjMapFrom1,)
TIMER(applyProjMapTo1,)
TIMER(applyProjMapFrom2,)
TIMER(applyProjMapTo2,)
void DerivativeFlat::ProjectionProduct::apply(){

    Parallel::debugTimer=false;
    ParallelProcess::debugTimer=false;
    for(size_t i=0; i<mapFrom.size(); i++){

        // These operations appear quite heavy (at least for a small system
        // some take equally long as the application of mapFrom/To).
        localSpec[i]->setToZero();
        localTemp->setToZero();
        globalSpec[i]->setToZero();
        globalTemp->setToZero();

        if(i==0) STARTDEBUG(applyProjMapFrom1) else STARTDEBUG(applyProjMapFrom2);
        mapFrom[i]->par->apply(1.);
        if(i==0) STOPDEBUG(applyProjMapFrom1) else STOPDEBUG(applyProjMapFrom2);
        if(i==0) STARTDEBUG(applyProjMapTo1) else STARTDEBUG(applyProjMapTo2);
        mapTo[i]->par->apply(1.);
        if(i==0) STOPDEBUG(applyProjMapTo1) else STOPDEBUG(applyProjMapTo2);

        // mapTo->par->apply(-1.) directly into cXY does not seem to work!
        *localXY -= *localTemp;
    }
    Parallel::debugTimer=false;
    ParallelProcess::debugTimer=false;

}

void DerivativeFlat::testProjection(){
    if(projection==0) return;

    MPIwrapper::Barrier();

    Coefficients c(setupXY->idx());
    c.setToRandom();

    Coefficients c1(setupXY->idx());
    Coefficients c2(setupXY->idx());

    Parallel::scatter(CoefficientsGlobal::view(&c), localXY, MPIwrapper::master());

    projection->apply();
    Parallel::gather(CoefficientsGlobal::view(&c1), localXY, MPIwrapper::master());

    projection->apply();
    Parallel::gather(CoefficientsGlobal::view(&c2), localXY, MPIwrapper::master());

    if(MPIwrapper::isMaster()){
        c1 -= c2;
        if(c1.norm()/c2.norm() > 1.e-9){
            std::cerr<<c1.str(0)<<std::endl;
            std::cerr<<"---------------------"<<std::endl;
            std::cerr<<c2.str(0)<<std::endl;
            std::cerr<<c1.norm()<<" / "<<c2.norm();
            DEVABORT("Failed to verify projector property");
        }
    }

}

TIMER(setup,)
void DerivativeFlat::_construct(const OperatorTree* Op, const DiscretizationSpectral *ProjectionDisc,
                                std::shared_ptr<ProjectSubspace> Project)
{
    STARTDEBUG(setup);
    timeCritical::suspend();

    _projSub=Project; //need to keep pointer alive...

    o=Op;
    iIndex=Op->iIndex;
    jIndex=Op->jIndex;

    inverseOverlap=Op->iIndex->inverseOverlap();

    setupTemp=0;
    globTemp=0;
    localTemp=0;

    setupXY=0;
    globXY=0;
    localXY=0;

    if(Op->iIndex!=Op->jIndex)ABORT("no flat derivative for operator between different spaces: "+Op->name);
    setupXY=  new Coefficients(Op->jIndex);
    setupTemp=new Coefficients(Op->jIndex);

    PrintOutput::message("running parallel with "+tools::str(MPIwrapper::Size())+" processes");

    LOG_PUSH("FlattenHamiltonian");
    op = new FlattenedOperatorTree(Op, true, applyEpsilon, setupXY, setupTemp, "either");

    if(op->par->str().find("(unused)")!=string::npos)PrintOutput::warning("unused process - choose fewer processes");
    LOG_POP();

    globTemp=CoefficientsGlobal::view(setupTemp);
    localTemp=CoefficientsLocal::view(globTemp);
    globXY=CoefficientsGlobal::view(setupXY);
    localXY=CoefficientsLocal::view(globXY);

    LOG_PUSH("FlattenProjection");
    projection=0;
    if(ProjectionDisc!=0){
        if(dynamic_cast<const DiscretizationSpectralProduct*>(ProjectionDisc) == 0){
            projection.reset(new ProjectionSingle(ProjectionDisc, setupXY, applyEpsilon));
        }else{
            projection.reset(new ProjectionProduct(dynamic_cast<const DiscretizationSpectralProduct*>(ProjectionDisc),
                                                   setupXY, applyEpsilon));
        }
    }
    else if(Project!=0){
        projection.reset(new ProjectionSingle(Project.get(), setupXY, applyEpsilon));
    }
    LOG_POP();

    LOG_PUSH("Prallel::test");

    Parallel::test(setupXY);
    LOG_POP();
    MPIwrapper::Barrier();

    // prepare for parallel application of inverse overlap
    globXY->idx()->inverseOverlap()->parallelSetup();

    // test projection
    LOG_PUSH("testProjection");
    testProjection();
    LOG_POP();
    timeCritical::resume();
    STOPDEBUG(setup);

}

void DerivativeFlat::project(Coefficients & C) const {
    if(not projection)return;

    *globXY=C;
    projection->apply();
    C=*globXY;
}


void DerivativeFlat::apply(std::complex<double> A, const Coefficients &X, std::complex<double> B, Coefficients &Y) const{
    //    ABORT("out of use");
    *globXY=X;
    apply(A,localXY,0.,*localXY);
    Y.scale(B);
    Parallel::gather(globTemp,localXY,0);
    if(MPIwrapper::isMaster())Y+=*globTemp;
}

TIMER(derPre,)
TIMER(derPre0,)
TIMER(derPre1,)
TIMER(derPre2,)
TIMER(derOp,)
TIMER(derInv,)
TIMER(derProj,)
TIMER(applyA1,)
TIMER(applyA2,)
TIMER(applyA3,)
TIMER(applyA4,)
TIMER(applyA5,)
TIMER(applyB1,)
TIMER(applyB2,)

void DerivativeFlat::applyA(std::complex<double> A, CoefficientsLocal *localX) const{
    //NOTE: global views should be replaced by halo-views

#ifdef _DEVELOP_
    if(localXY->isNan()){
        Mout+localX->str()+Sendl;
        DEVABORT("local nan A in");
    }
#endif
    STARTDEBUG(applyA1);
    op->par->setToZeroLHS();
    STOPDEBUG(applyA1);

    STARTDEBUG(applyA2);
    *localTemp=*localX;
    STOPDEBUG(applyA2);

    STARTDEBUG(applyA3);

    op->par->apply(A*complex<double>(0.,-1.));
    STOPDEBUG(applyA3);

    STARTDEBUG(applyA4);
    *localTemp=*localXY;
    STOPDEBUG(applyA4);

    STARTDEBUG(applyA5);
    inverseOverlap->apply(1.,*localTemp,0.,*localXY);
    STOPDEBUG(applyA5);
#ifdef _DEVELOP_
    if(localXY->isNan()){
        Mout+"temp"+localTemp->isNan()+Sendl;
        Mout+"local"+localXY->isNan()+Sendl;
        MPIwrapper::Barrier();
        DEVABORT("local nan A out");
    }
#endif
}
void DerivativeFlat::applyB(std::complex<double> B, CoefficientsLocal &Y) const{

    // CAUTION:
    // for efficiency reasons, initial state is not projected
    // if initial state is not in proper subspace,
    // non-hermitian operators and norm-nonconservation results
    STARTDEBUG(applyB1);
    if(projection!=0)projection->apply();
    STOPDEBUG(applyB1);

    STARTDEBUG(applyB2);
    if(B==0.)
        Y=*localXY;
    else {
        if(B!=1.)Y*=B;
        Y+=*localXY;
    }
    STOPDEBUG(applyB2);

}
void DerivativeFlat::apply(std::complex<double> A, CoefficientsLocal *localX, std::complex<double> B, CoefficientsLocal &Y) const{
    applyA(A,localX);
    applyB(B,Y);
}

void DerivativeFlat::eigenValues(double Time, std::vector<complex<double> > Eval) {
    ABORT("need to fix - apply overlap");
    UseMatrix mat;
    matrix(mat);
    o->iIndex->matrixContract(mat,mat);
    mat*=complex<double>(0.,1.);
    vector<complex<double> > ovr(idx()->dvrWeigContract());
    for(unsigned int k=0;k<mat.cols();k++)mat.col(k)/=sqrt(ovr[k]);
    for(unsigned int k=0;k<mat.rows();k++)mat.row(k)/=sqrt(ovr[k]);
    UseMatrix eval;
    mat.eigenValues(eval);
    Eval.clear();
    vector<vector<double> > realImag(2,vector<double>());
    for(unsigned int k=0;k<eval.size();k++){
        Eval.push_back(eval(k).complex());
        realImag[0].push_back(Eval.back().real());
        realImag[1].push_back(Eval.back().imag());
    }
    AsciiFile eig(ReadInput::main.output()+"flatEigen");
    eig.writeCols(realImag);
    PrintOutput::message("eigenvalues on "+eig.name());
}

DerivativeFlat::Arp::Arp(const DerivativeFlat *Der)
    : Arpack(1.e-9,5000),der(Der),x(Der->lhsVector().idx()),y(Der->idx()),
      locX(CoefficientsLocal(Der->idx())),locY(CoefficientsLocal(Der->idx()))
{
    x.pointerToC(pX);
    y.pointerToC(pY);
    _lvec=x.size();
}

void DerivativeFlat::eigen(std::vector<std::complex<double> > & Eval,std::vector<Coefficients*> &Evec, unsigned int NStat){
    //NOTE: this should be funneled through an EigenSolverAbstract
    Arp A(this);
    string method="SmallReal";
    if(Evec.size()==1)method="SmallAbs";
    A.eigen(Eval,Evec,NStat,method,Evec.size()==1);

    // fix the phases for reproducible results
    for(Coefficients* c: Evec){
        std::complex<double> phase=c->cMaxNorm();
        c->scale(std::conj(phase/std::abs(phase)));
    }
}

static int __counter = 0;

void DerivativeFlat::Arp::eigen(vector<std::complex<double> > &Eval, vector<Coefficients*> &Rvec,
                                unsigned int Nvec, const string &Which, bool Restart){

    LOG_MEM("Before DerivativeFlat::Arp::eigen");
    LOG_PUSH("DerivativeFlat::Arp::eigen");
    __counter = 0;

    complex<double> defaultShift=shift;
    if(Restart){
        // initialize with starting vector
        if(Rvec.size()!=1)ABORT("for restarting, supply exactly one Rvec as starting vector");
        rvec.resize(_lvec);
        for(unsigned int i=0;i<_lvec;i++)rvec[i]=Rvec[0]->data()[i];
        delete Rvec[0];
        Rvec.clear();
        shift=Eval[0];
    }
    else if(Rvec.size()>0)
        ABORT("must enter eigen with empty Rvec, unless Restart");

    if(Which=="SmallAbs"){
        if(Eval.size()!=1)ABORT("for Which=SmallAbs provide guess energy value in Eval");
        shift=Eval[0];
    }

    // let ArpackFunction compute the vectors
    Arpack::eigenIter(Nvec,Which,Restart,true);

    // move into Rvec
    Eval.resize(Nvec);
    for(unsigned int k=0;k<Nvec;k++){
        Eval[k]=eval[k]+shift;
        Rvec.push_back(new Coefficients(y.idx()));
        for(size_t i=0;i<Rvec.back()->size();i++)Rvec.back()->data()[i]=rvec[i+k*_lvec];
        *Rvec.back()*=1./sqrt(der->idx()->overlap()->matrixElement(*Rvec.back(),*Rvec.back(),false));
    }

    // reset in case is had been changed
    shift=defaultShift;
    LOG_POP();
}

void DerivativeFlat::Arp::apply(const std::complex<double> *X, std::complex<double> *Y) {
    __counter++;
    if(__counter%10 == 0) LOG_D(__counter);

    // copy to Coefficients
    for (unsigned int i=0;i<pX.size();i++)*(pX[i])=*(X+i);
    if(MPIwrapper::isMaster())x.makeContinuous();
    // get zero-projected content and shift
    for (unsigned int i=0;i<pX.size();i++)Y[i]=shiftZero*(X[i]-*(pX[i]));

    Parallel::scatter(&x,&locX,MPIwrapper::master());

    // y=Der(x)-shift*x
    der->apply(complex<double>(0.,1.),&locX,0.,locY);
    locY.axpy(-shift,&locX);
    // copy back and reset y to zero
    Parallel::gather(&y,&locY,MPIwrapper::master());
    for (unsigned int i=0;i<pX.size();i++){
        *(Y+i)+=*(pY[i]); // contains shifted zero-projected content
        *(pY[i])=0.;
    }
}
