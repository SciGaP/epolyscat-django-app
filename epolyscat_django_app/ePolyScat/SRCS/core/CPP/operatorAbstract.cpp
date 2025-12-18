// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "operatorAbstract.h"

#include <stdlib.h>

#include "printOutput.h"
#include "plot.h"
#include "index.h"
#include "inverse.h"
#include "coefficients.h"
//#include "operator.h"
#include "operatorTree.h"
#include "operatorMap.h"
#include "mapGridHybrid.h"

#include "overlapDVR.h"
#include "constrainedView.h"

#include "timer.h"
#include "str.h"
#include "parameters.h"
#include "timeCritical.h"
#include "readInput.h"
#include "operatorDefinition.h"

#include "tools.h"
#include "folder.h"
#include "blockView.h"

using namespace std;

bool OperatorAbstract::useOperatorFloor=true;
bool OperatorAbstract::fuseOp=true;
bool OperatorAbstract::useTensor=false;
bool OperatorAbstract::flat=true;
string OperatorAbstract::eigenMethod="auto";          ///< allow arnoldi eigensolver (not stable at times)
std::complex<double> OperatorAbstract::Arp_shift=complex<double>(100.);
void OperatorAbstract::readControls(ReadInput &Inp){
    EigenSolverAbstract::readControls(Inp);
    OperatorAbstract::eigenMethod=EigenSolverAbstract::defaultMethod;
    double inShift;
    Inp.obsolete("Eigensolver","shift","use Eigen: shift instead");
    Inp.read("Eigen","shift",inShift,"100","shift Arpack eigensolver to avoid zero's arising from continuity"
                                           " projections if positive eigenvalues are searched for")
            .texdocu(R"exp(
                     When iterative methods are use, continuity constraints lead to zero eigenvalues in the full Hamiltonian
                     (which correspond to discontinous solutions). When looking for a positive eigenalues, this shifts
                     the hamiltonian to negative values (while the continuity zeros remain at 0), which will be found as the lowest
                     eigenvalues of the problem. The shift will be automatically deducted from the final result.
                     With a default of 100 there is usually no need to specify this.
                     )exp");
    OperatorAbstract::Arp_shift=inShift;
    ReadInput::main.read("_EXPERT_","newFloor",useOperatorFloor,"false","admit use of tensor products",0,"newFloor");
    ReadInput::main.read("_EXPERT_","suppressTensor",OperatorAbstract::useTensor,"false","admit use of tensor products",0,"noTensor");
    OperatorAbstract::useTensor=not OperatorAbstract::useTensor;

    Inp.read("_EXPERT_","fuseOperator",fuseOp,"true","fuse pieces of the operator as much as possible",1,"fuse");
    OperatorAbstract::useTensor=useTensor and (not fuseOp);

    Inp.read("_EXPERT_","flatDer",OperatorAbstract::flat,"true","absorb the inverse into new Floor operators",1,"flatDer");
    if(flat){
        //        OperatorFloor::absorbInverse=true;
        useOperatorFloor=true;
        useTensor=false;
    }

    PrintOutput::DEVmessage("Expert controls: fuse="+tools::str(fuseOp)
                            +", useTensor="+tools::str(useTensor)
                            +", newFloor="+tools::str(useOperatorFloor)
                            +", flatDerivative="+tools::str(OperatorAbstract::flat)
                            //                         +", absorbInverse="+tools::str(OperatorFloor::absorbInverse)
                            );
}

static void writeString(std::ofstream *File,const std::string S){
    tools::write(*File,int(S.length()));
    tools::write(*File,S.data(),S.length());
}
static std::string readString(std::ifstream &File){
    int siz;
    tools::read(File,siz);
    char s[siz];
    tools::read(File,s,siz);
    return std::string(s,siz);
}

Index* OperatorAbstract::idx(bool True){
    return const_cast<Index*>(True?iIndex:jIndex);
}

//std::string OperatorAbstract::plotActionOnCoefficients(const Coefficients * C, std::shared_ptr<Plot> Plt) const
//{
//    if(Plt==0)Plt.reset(new Plot(idx(),ReadInput::main));
//    if(Plt->isEmpty())return "cannot plot action of "+name+", neet Plot: ... in input";

//    Coefficients c(jdx()),d(idx()),p(idx());
//    if(C!=0)c=*C;
//    else c.setToFunction("1");

//    apply(1.,c,0.,d);
//    idx()->inverseOverlap()->apply(1,d,0.,p);

//    std::string file=name+(C==0?"_1":"_C");
//    file=ReadInput::main.output()+file;

//    Plt->withPlotReal().plot(p,file);
//    return "action of "+name+" on "+file;
//}

ofstream *OperatorAbstract::write(ofstream *File){
    if(not File->is_open())DEVABORT("cannot write - File not open");
    int code(5);
    tools::write(*File,code); // version code
    writeString(File,ReadInput::main.outputTopDir());
    writeString(File,name);
    writeString(File,definition.str());
    iIndex->write(*File,true);
    jIndex->write(*File,true);
    return File;
}
static void checkAndAssign(const Index* &Required, const Index* Read){
    if(Required){
        if(not Required->treeEquivalent(Read))
            ABORT("cached iIndex\n"+Read->str()+"\ndiffers from required\n"+Required->str());
        delete Read;
    }
    else Required=Read;
}
OperatorAbstract::OperatorAbstract(std::ifstream &File, string Name, const Index *IIndex, const Index* JIndex)
    :OperatorAbstract(Name,IIndex,JIndex){
    if(not File.is_open())ABORT("input stream not open - check file name");
    if(File.tellg()==0){
        int code;
        tools::read(File,code);
        if(code==5){
            std::string created_by=readString(File);
            name=readString(File);
            if(Name!="" and name!=Name)ABORT("cached name "+name+" differs from required "+Name);
            definition=OperatorDefinition(readString(File));
            checkAndAssign(iIndex,new Index(File,false));
            const Index* jdx=new Index(File,false);
            if(jdx->treeEquivalent(iIndex))jIndex=iIndex;
            else checkAndAssign(jIndex,jdx);
        }
        else
            DEVABORT(Sstr+"unknown operator storage code: "+code);

    }
}

OperatorAbstract::~OperatorAbstract(){
    delete _tempRHS;
    delete _tempLHS;
}

Coefficients* OperatorAbstract::tempLHS() const {
    if(iIndex==0)return 0;
    if(_tempLHS==0){
        _tempLHS=new Coefficients(iIndex);
    }
    return _tempLHS;
}
Coefficients* OperatorAbstract::tempRHS() const {
    if(jIndex==0)return 0;
    if(_tempRHS==0){
        timeCritical::suspend(); // exempted from time-critical warnings
        _tempRHS=new Coefficients(jIndex);
        timeCritical::resume();
    }
    return _tempRHS;
}

void OperatorAbstract::axpy(std::complex<double> Alfa,const Coefficients& X,std::complex<double> Beta, Coefficients& Y, const double Time)
{
    update(Time);
    apply(Alfa,X,Beta,Y);
}

bool OperatorAbstract::isIdentity(double Eps, bool Stochastic) const {
    tempRHS()->setToRandom();
    *tempLHS()=*tempRHS();
    apply(-1.,*tempRHS(),1.,*tempLHS());
    if(not tempLHS()->isZero(Eps))return false;
    if(not Stochastic) DEVABORT("only stochastic test implemented");
    return true;
}

bool OperatorAbstract::isZero(double Eps, bool Stochastic) const {
    tempRHS()->setToRandom();
    apply(-1.,*tempRHS(),0.,*tempLHS());
    if(not tempLHS()->isZero(Eps))return false;
    if(not Stochastic) DEVABORT("only stochastic test implemented");
    return true;
}

bool OperatorAbstract::isSelfAdjoint(double Eps) const{return isSymmetric("SelfAdjoint",Eps);}
bool OperatorAbstract::isComplexSymmetric(double Eps) const{return isSymmetric("ComplexSymmetric",Eps);}
bool OperatorAbstract::isSymmetric(std::string Kind,double Eps) const{
    const OperatorTree* oTree=dynamic_cast<const OperatorTree*>(this);
    bool res=true;
    if(iIndex->size()!=jIndex->size())
        res=false;
    else if(oTree)
        res=oTree->isSymmetric(Kind,Eps);
    else {
        if(isHuge()){
            Coefficients l(idx()),r(jdx());
            l.setToRandom();
            r.setToRandom();
            std::complex<double> lr=matrixElement(l,r);
            if(Kind=="ComplexSymmetric")
                res=std::abs(lr-matrixElement(r,l))>std::abs(lr)*1.e-10;
            else if(Kind=="SelfAdjoint"){
                l.conjugate();
                r.conjugate();
                res=std::abs(std::conj(lr)-matrixElement(r,l))>std::abs(lr)*1.e-10;
            }
            else DEVABORT("illegal Kind="+Kind);
        }
        else {
            Eigen::MatrixXcd mat=matrix();
            for(int k=0;k<mat.rows();k++){
                if(Kind=="SelfAdjoint"){
                    if(!((mat.col(k)-mat.row(k).adjoint()).isZero(Eps*mat.col(k).lpNorm<Eigen::Infinity>())))return false;
                }
                else if(Kind=="ComplexSymmetric"){
                    if(!(mat.col(k)-mat.row(k).transpose()).isZero(Eps*mat.col(k).lpNorm<Eigen::Infinity>()))return false;
                }
                else
                    DEVABORT("illegal Kind="+Kind);
            }
        }
    }
    return res;
}


UseMatrix OperatorAbstract::matrix(std::vector<Coefficients *> Left, std::vector<Coefficients *> Right) const{
    Coefficients opTimesRight(jIndex);
    UseMatrix res(Left.size(),Right.size());
    for(size_t j=0;j<Right.size();j++){
        apply(1.,*Right[j],0.,opTimesRight);
        for(size_t i=0;i<Left.size();i++)
            res(i,j)=Left[i]->innerProduct(&opTimesRight,false);
    }
    return res;
}
Eigen::MatrixXcd OperatorAbstract::matrix() const{
    if(isHuge())ABORT("too big - cannot construct full matrix");
    vector<complex<double>> mat;
    matrix(mat);
    return Eigen::Map<Eigen::MatrixXcd>(mat.data(),iIndex->size(),jIndex->size());
}

Eigen::SparseMatrix<std::complex<double>> OperatorAbstract::matrixSparse(bool Contract, double Eps) const {
    if(isHuge())ABORT("huge matrix");
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
    }

    // convert to triples wrt to global index, normalize where multiplicities
    idim=0;
    jdim=0;
    Parameters::updateToOne(); // set all time-dependent parameters =1
    //    for(size_t k=0;k<_trip.size();k++){
    Coefficients rhs(jdx());
    Coefficients lhs(idx());
    std::vector<double> maxCol(jdim);
    for(size_t j=0;j<jdx()->size();j++){
        rhs.data()[j]=1.;
        apply(1.,rhs,0.,lhs);
        maxCol[jglob[j]]=std::max(maxCol[jglob[j]],(Eigen::Map<Eigen::MatrixXcd>(lhs.data(),lhs.size(),1).lpNorm<Eigen::Infinity>()));
        for(size_t i=0;i<lhs.size();i++){
            complex<double> matij=lhs.data()[i];
            if(inrm.size()>0)matij*=inrm[i];
            if(jnrm.size()>0)matij*=jnrm[j];
            idim=max(idim,(int)iglob[i]);
            jdim=max(jdim,(int)jglob[j]);
            if(std::norm(matij)>0.)list.push_back(Eigen::Triplet< complex<double> >(iglob[i],jglob[j],matij));
        }
        rhs.data()[j]=0.;

    }
    Parameters::restoreToTime();
    idim++;
    jdim++;

    // remove near-zeros from triplets
    // criterion uses only column-norms
    // for symmetric/hermitian matrices criterion matches block-view criterion
    double maxAll=*std::max(maxCol.begin(),maxCol.end());
    size_t nzero=0;
    for(size_t k=0;k<list.size();k++){
        if(std::norm(list[k].value())>
                Eps*Eps*(list[k].row()<(jdim?maxCol[list[k].row()]:maxAll)*maxCol[list[k].col()])){
            list[nzero]=list[k];
            nzero++;
        }
    }
    list.resize(nzero);
    Eigen::SparseMatrix< complex<double>> res(idim,jdim);
    res.setFromTriplets(list.begin(),list.end(),[] (const  complex<double> & a,const  complex<double> &b) { return a+b; });

    return res;
}

bool OperatorAbstract::isHuge(std::string Message) const{
    if(ReadInput::main.flag("DEBUGallowHugeMatrix","suppress checks on huge matrices"))return false;
    if(iIndex->size()*jIndex->size()>4000*4000){
        PrintOutput::warning(Sstr+"huge matrix"+name+"dimensions"+iIndex->size()+"x"+jIndex->size()
                             +"suppress check by command line option -DEBUGallowHugeMatrix\n"+Message);
        return true;
    }
    return false;
}

void OperatorAbstract::matrix(UseMatrix &Mat) const{
    if(isHuge())ABORT("too big - cannot construct full matrix");
    vector<complex<double> > mat;
    matrix(mat);
    Mat=UseMatrix::UseMap(mat.data(),iIndex->sizeCompute(),jIndex->sizeCompute());
}

TIMER(matrix,)
void OperatorAbstract::matrix(std::vector<std::complex<double> > &Mat) const{
    if(isHuge())ABORT("too big - cannot construct full matrix");
    STARTDEBUG(matrix);
    if(Parameters::currentTime()==-DBL_MAX){
        Parameters::updateToOne();
        const OperatorTree * oTree=dynamic_cast<const OperatorTree *>(this);
        if(oTree==0 or oTree->def().find("[t]")!=string::npos)
            PrintOutput::DEVwarning("no time has been set - setting all time-dependent parameters =1");
    }

    iIndex->sizeCompute();
    jIndex->sizeCompute();
    Coefficients a(iIndex),b(jIndex);
    a.treeOrderStorage();
    b.treeOrderStorage();
    Mat.assign(iIndex->sizeStored()*jIndex->sizeStored(),0.);
    for(size_t j=0,ij=0;j<jIndex->sizeStored();j++){
        b.setToZero();
        a.setToZero();
        b.storageData()[j]=1.;

        apply(1.,b,0.,a);
        for(size_t k=0;k<iIndex->sizeStored();k++,ij++){
            Mat[ij]=a.storageData()[k];
        }
    }
    Parameters::restoreToTime();

    STOPDEBUG(matrix)
}

Eigen::MatrixXcd OperatorAbstract::matrixContracted() const{
    vector<complex<double> > mat;
    matrix(mat);
    UseMatrix Mat=UseMatrix::UseMap(mat.data(),iIndex->sizeStored(),jIndex->sizeStored());
    UseMatrix res;
    matrixContract(Mat,res,1.);
    return Eigen::Map<Eigen::MatrixXcd>(res.data(),res.rows(),res.cols());
}

void OperatorAbstract::matrixAdd(std::complex<double> factor, UseMatrix &Mat) const{
    vector<complex<double> > mat;
    matrix(mat);
    Mat=UseMatrix::UseMap(mat.data(),iIndex->sizeStored(),jIndex->sizeStored());
    matrixContract(Mat,Mat,factor);
}

complex<double> OperatorAbstract::matrixElement(const Coefficients &Ci, const Coefficients &Cj, bool pseudoScalar) const{
    apply(1.,Cj,0.,*tempLHS());
    tempLHS()->makeContinuous();
    return Ci.innerProduct(tempLHS(),pseudoScalar);
}
complex<double> OperatorAbstract::matrixElementUnscaled(const Coefficients &Ci, const Coefficients &Cj) const{
    apply(1.,Cj,0.,*tempLHS());
    complex<double> tmp=Ci.innerProductUnscaled(tempLHS());
    return tmp;
}

TIMER(contCol,)
TIMER(contRow,)
TIMER(contMult,)

/// GMat with continuity imposed (Mat and GMat may be the same UseMatrix)
void OperatorAbstract::matrixContract(const UseMatrix &Mat,UseMatrix & GMat,complex<double> Factor,std::vector<int>ISort,std::vector<int>JSort) const {
    // get global indices
    vector<unsigned int>iglob,jglob;
    std::vector<double> inorm,jnorm;
    idx()->multiplicities(iglob,inorm);
    jdx()->multiplicities(jglob,jnorm);

    // check matrix dimensions
    if(iglob.size()!=Mat.rows() or jglob.size()!=Mat.cols())
        ABORT("input matrix dimensions do not match index size "
              +tools::str(int(iglob.size()))+" X "+tools::str(int(jglob.size()))+" vs. "+tools::str(Mat.rows())+" "+tools::str(Mat.cols()));

    STARTDEBUG(contCol);
    // contract columns
    UseMatrix cmat=UseMatrix::Zero(Mat.rows(),jIndex->globalLength());
    if(JSort.size()>0 and JSort.size()!=cmat.cols())ABORT("JSort does not match global length");
    for (unsigned int j=0;j<Mat.cols();j++){
        unsigned int js=jglob[j];
        if(JSort.size()>0)js=JSort[jglob[j]];
        cmat.block(0,js,Mat.rows(),1)+=Mat.col(j)*jnorm[j];
    }
    STOPDEBUG(contCol);

    STARTDEBUG(contRow);
    // contract rows
    unsigned int isiz=*std::max_element(iglob.begin(),iglob.end())-*std::min_element(iglob.begin(),iglob.end())+1;
    unsigned int jsiz=*std::max_element(jglob.begin(),jglob.end())-*std::min_element(jglob.begin(),jglob.end())+1;
    GMat=UseMatrix::Zero(isiz,jsiz);
    if(ISort.size()>0 and ISort.size()!=GMat.rows())ABORT("ISort does not match global length");
    for (unsigned int i=0;i<cmat.rows();i++){
        unsigned int is=iglob[i];
        if(ISort.size()>0)is=ISort[iglob[i]];
        for(size_t k=0,kg=is,kc=i;k<cmat.cols();k++,kg+=GMat.leadDim(),kc+=cmat.leadDim())
            GMat.data()[kg]+=cmat.data()[kc]*Factor*inorm[i];
    }
    STOPDEBUG(contRow);
}

bool OperatorAbstract::isBlockDiagonal() const {
    if(dynamic_cast<const ConstrainedView*>(this)!=0)dynamic_cast<const ConstrainedView*>(this)->isBlockDiagonal();
    if(dynamic_cast<const OverlapDVR*>(this)!=0)return iIndex->continuity()==Index::npos;
    if(dynamic_cast<const OperatorTree*>(this)!=0)return dynamic_cast<const OperatorTree*>(this)->isBlockDiagonal();
//    if(dynamic_cast<const Operator*>(this)!=0)return dynamic_cast<const Operator*>(this)->isBlockDiagonal();
    return false;
}

bool OperatorAbstract::isDiagonal() const {
    DEVABORT("no generic implemementation");
    return false;
}

void OperatorAbstract::subMatrix(std::vector<std::complex<double> > &Mat, const Index *ISub, const Index *JSub) const{
    Coefficients a(iIndex),b(jIndex);
    a.treeOrderStorage();
    b.treeOrderStorage();
    Mat.resize(ISub->sizeStored()*JSub->sizeStored());

    for(size_t j=JSub->posIndex(jIndex),ij=0;j<JSub->posIndex(jIndex)+JSub->sizeStored();j++){
        b.setToZero();
        b.data()[j]=1.;
        apply(1.,b,0.,a);
        for(size_t i=ISub->posIndex(iIndex);i<ISub->posIndex(iIndex)+ISub->sizeStored();i++,ij++){
            Mat[ij]=a.data()[i];
        }
    }
}

void OperatorAbstract::orthoNormalizeDegenerate(std::vector<std::complex<double> > Eval, std::vector<Coefficients *> Evec, bool Pseudo) const{
    double epsRelDegen=1.e-12,epsAbsDegen=1.e-8; // level where to consider eigenvalues degenerated

    if(Eval.size()!=Evec.size())ABORT("number of eigenvalues and eigenvectors does not match");

    // orthonormalize (near-)degenerate subspaces
    for(size_t k=0;k<Evec.size();k++){
        vector<Coefficients*>sub(1,Evec[k]);
        for(size_t l=k+1;l<Eval.size();l++){
            if(abs(Eval[l]-Eval[k])<epsRelDegen*abs(Eval[k])
                    or abs(Eval[l]-Eval[k])<epsAbsDegen)
                sub.push_back(Evec[l]);
        }
        orthoNormalize(sub,Pseudo);
    }
}

void OperatorAbstract::orthoNormalize(std::vector<Coefficients *> & C, bool pseudo) const {
    if(C.size()==0)return;
    Coefficients *sCk=tempLHS();
    for(size_t k=0;k<C.size();k++){
        apply(1.,*C[k],0.,*sCk);
        for(int l=0;l<k;l++)C[k]->axpy(-C[l]->innerProduct(sCk,pseudo),C[l]);
        complex<double> a=C[k]->innerProduct(sCk,pseudo);
        if(abs(a)<1.e-12){
            if(abs(a)==0.)ABORT("vectors linearly dependent");
            PrintOutput::warning("near-linearly dependent vectors");
        }
        *C[k]*=1./sqrt(a);
    }
}

// TODO: Copied from operatorFloor.cpp...
static Timer TIMERcost("cost","INTERNAL","operatorAbstract.cpp");
double OperatorAbstract::applicationCost() const{

    Coefficients x(jIndex);
    Coefficients y(iIndex);

    x.treeOrderStorage();
    y.treeOrderStorage();

    x.setToRandom();
    y.setToRandom();

    TIMERcost.start();
    TIMERcost.stopTimer();
    double startSecs=TIMERcost.secs();
    TIMERcost.start();

    for(unsigned int k=0;k<10;k++){
        apply(1.,x,0.,y);
        apply(1.,x,0.,y);
        apply(1.,x,0.,y);
        apply(1.,x,0.,y);
        apply(1.,x,0.,y);
        apply(1.,x,0.,y);
        apply(1.,x,0.,y);
        apply(1.,x,0.,y);
        apply(1.,x,0.,y);
        apply(1.,x,0.,y);
    }

    TIMERcost.stopTimer();
    return 1.e5*(TIMERcost.secs()-startSecs);

}

double OperatorAbstract::norm() const{ABORT("not implemented");}

std::string OperatorAbstract::str(int Level) const
{
    if(dynamic_cast<const OperatorTree*>(this))
        return dynamic_cast<const OperatorTree*>(this)->str(Level);
    else if(dynamic_cast<const OperatorMap*>(this))
        return dynamic_cast<const OperatorMap*>(this)->Tree::str(Level);
    else if(dynamic_cast<const MapGridHybrid*>(this))
        return dynamic_cast<const MapGridHybrid*>(this)->Tree::str(Level);
    else {
        Str s("");
        s=s+iIndex->index()+" | "+jIndex->index();
        if(Level==Tree_defaultKind)return s+name+": "+iIndex->hierarchy()+" <-- "+jIndex->hierarchy();
        else return s+Level;
    }
}
