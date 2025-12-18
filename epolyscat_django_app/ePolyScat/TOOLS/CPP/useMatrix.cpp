// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "useMatrix.h"
#include "tools.h"
#ifdef _HAS_LAPACKE_
#include "lapacke.h"
#endif
#ifdef _USE_ZSLAPACK_
#include "zslapack.h"
#endif
#include "abort.h"
#include "triFactorLU.h"
#include "printOutput.h"
#include <Eigenvalues>
#include "qtEigenDense.h"

#include "timer.h"
TIMER(eigen,)
TIMER(compress,)
TIMER(solve,)

using namespace std;

// convenient short hand notation
typedef complex<double>* Z;
typedef complex<double> z;
typedef double* R;
typedef double r;

static const unsigned int smallDimensions=20;

// access to data
// basic storage formats
unsigned int UseMatrix::full_normal(unsigned int I,unsigned int J, const Shape &shape){return I+J*shape.leadDim;}
unsigned int UseMatrix::band_normal(unsigned int I,unsigned int J, const Shape &shape){
    if(I<J){if(J-I+1>shape.upBw )return shape.total-1;}
    else {  if(I-J+1>shape.lowBw)return shape.total-1;}
    return shape.upBw-1+I-J+J*shape.leadDim;
}
unsigned int UseMatrix::band_symmetric(unsigned int I,unsigned int J, const Shape & shape){
    if(I>J){
        if(I-J > shape.lowBw-1)return shape.total-1;
        return I-J+J*shape.leadDim;
    } else {
        if(J-I > shape.upBw -1)return shape.total-1;
        return J-I+I*shape.leadDim;
    }
}
// views:
// transpose
unsigned int UseMatrix::any_transpose  (unsigned int I,unsigned int J, const Shape &shape)
{
    //    cout<<"transposed J,I "<<J<<" "<<I<<": "<< shape.stored->locationIJ(J,I,*(shape.stored))<<"..."<<shape.stored->str()<<endl;
    return shape.stored->locationIJ(J,I,*(shape.stored));}
// block
unsigned int UseMatrix::any_block   (unsigned int I,unsigned int J, const Shape &shape)
{return shape.stored->locationIJ(shape.row0+I,shape.col0+J,*(shape.stored));}
// band
unsigned int UseMatrix::any_band   (unsigned int I,unsigned int J, const Shape &shape)
{   if(I-J+1>shape.lowBw or J-I+1>shape.upBw)return shape.stored->total-1; // outside band - point to last (=0) element
    return shape.stored->locationIJ(I,J,*(shape.stored));}
// map full
unsigned int UseMatrix::full_normal_map   (unsigned int I,unsigned int J, const Shape &shape)
{return full_normal(I,J,shape);}
// map banded
unsigned int UseMatrix::asDiagonal_map   (unsigned int I,unsigned int J, const Shape &shape)
{return band_normal(I,J,shape);}

// in place scalar operations, e.g. A*=a;
const UseMatrix::Operation  UseMatrix::ScalarInAssign =UseMatrix::Operation( assign,scalar,true);
const UseMatrix::Operation  UseMatrix::ScalarInAdd    =UseMatrix::Operation(    add,scalar,true);
const UseMatrix::Operation  UseMatrix::ScalarInProduct=UseMatrix::Operation(product,scalar,true);

// out-of-place scalar operations, e.g. return  A+a;
const UseMatrix::Operation  UseMatrix::ScalarOutAdd    =UseMatrix::Operation(    add,scalar,false);
const UseMatrix::Operation  UseMatrix::ScalarOutProduct=UseMatrix::Operation(product,scalar,false);

// in place element-wise binary operations, e.g. A+=B;
const UseMatrix::Operation UseMatrix::BinaryInAssign  =UseMatrix::Operation(  assign,binary,true);
const UseMatrix::Operation UseMatrix::BinaryInAdd     =UseMatrix::Operation(     add,binary,true);
const UseMatrix::Operation UseMatrix::BinaryInSubtract=UseMatrix::Operation(subtract,binary,true);
const UseMatrix::Operation UseMatrix::BinaryInProduct =UseMatrix::Operation( product,binary,true);
const UseMatrix::Operation UseMatrix::BinaryInQuotient=UseMatrix::Operation(quotient,binary,true);

// out-of-place ternary operations, e.g. C=A+B
const UseMatrix::Operation UseMatrix::Add     =UseMatrix::Operation(     add,ternary,false);
const UseMatrix::Operation UseMatrix::Subtract=UseMatrix::Operation(subtract,ternary,false);
const UseMatrix::Operation UseMatrix::Product =UseMatrix::Operation( product,ternary,false);
const UseMatrix::Operation UseMatrix::Quotient=UseMatrix::Operation(quotient,ternary,false);

bool UseMatrix::Shape::operator==(const UseMatrix::Shape & B) const
{
    if(nrows!=B.nrows)return false;
    if(ncols!=B.ncols)return false;
    if(total!=B.total)return false;
    if(locationIJ!=B.locationIJ)return false;
    if(leadDim!=B.leadDim)return false;
    if(lowBw!=B.lowBw)return false;
    if( upBw!=B.upBw)return false;
    if((stored==0)!=(B.stored==0))return false;
    if(stored!=0 and !(stored->operator==(*B.stored)))return false;
    return true;
}


UseMatrix::~UseMatrix(){freeData();delete trf;delete symm;}
void UseMatrix::freeData(){
    // only standard matrices own their data
    if(isBasic()){
        delete[] cdata;
        delete[] rdata;
    }
    rdata=0;
    cdata=0;
    delete trf;
    trf=0;
}

string UseMatrix::Shape::strLocIJ() const {
    if(locationIJ==full_normal)return"full_normal";
    if(locationIJ==band_normal)return"band_normal";
    if(locationIJ==band_symmetric)return"band_symmetric";
    if(locationIJ==full_normal_map)return"full_normal_map";
    // views
    if(locationIJ==any_transpose)return "transpose("+stored->strLocIJ()+")";
    if(locationIJ==any_block)return "block("+stored->strLocIJ()+")";
    if(locationIJ==asDiagonal_map)return "diagonal("+stored->strLocIJ()+")";
    return "unknown conversion to string";
}

bool UseMatrix::isBasic() const{
    if(conjugated)return false;
    return shape.storeBasic();
}
bool UseMatrix::Shape::storeBasic() const {
    if(locationIJ==full_normal)return true;
    if(locationIJ==band_normal)return true;
    if(locationIJ==band_symmetric)return true;
    if(locationIJ==full_normal_map)return false;
    if(locationIJ==asDiagonal_map)return false;
    return stored==0; // matrices with storage different from shape are never basic
    ABORT("cannot decide basic: "+strLocIJ());
    abort();
}
bool UseMatrix::Shape::storeFull() const {
    if(locationIJ==band_normal       )return false;
    if(locationIJ==band_symmetric    )return false;
    if(stored!=0)return stored->storeFull();
    return true;
}
bool UseMatrix::Shape::storeTranspose() const {
    if(locationIJ==full_normal)return false;
    if(locationIJ==band_normal)return false;
    if(locationIJ==band_symmetric)return false;

    // views
    if(locationIJ==full_normal_map)return false;
    if(locationIJ==any_block)return stored->storeTranspose();
    if(locationIJ==any_transpose)return not stored->storeTranspose();
    if(locationIJ==asDiagonal_map)return false;
    cout<<locationIJ<<endl;
    ABORT("undefined transposition state: "+strLocIJ());
    abort();
}
bool UseMatrix::Shape::storeSymmetric() const {
    if(locationIJ==band_symmetric)return true;
    return false;
}

bool UseMatrix::Shape::storeContiguous() const{
    if(storeBasic())return true;
    if(locationIJ==any_transpose)return stored->storeContiguous();
    return false;
}

void UseMatrix::Shape::print(string text) const {
    if(text.length()>0)cout<<text<<": ";
    cout<<str()<<endl;
}

string UseMatrix::Shape::str() const {
    return strLocIJ()+": "+tools::str(nrows)+" x "+tools::str(ncols)+" (total,leadDim,row0,col0,lowBw,upBw= "
            +tools::str(total)+", "+tools::str(leadDim)+", "+tools::str(row0)+", "+tools::str(col0)
            +", "+tools::str(lowBw)+", "+tools::str(upBw)+")";
}

template<class D>
void UseMatrix::cornerData(const D & data) const {
    if(isFull())return;
    // unused corners of banded storage are sometimes accessed for speed, make sure they have legal values

    //lower right corner
    D end=data+shape.leadDim;
    for(unsigned int j=cols()-shape.lowBw+1;j<cols();j++,end+=shape.leadDim)
        for(D a=endCol(data,j);a<end;a++)*a=1.;

    if(isShape(band_symmetric))return; // upper left corner is not there

    // upper left corner
    D start=data;
    for(unsigned int j=0; j<shape.upBw-1;start+=shape.leadDim,j++)
        for(D a=start;a<begCol(data,j);a++)*a=1.;
}
// allocate data as needed
void UseMatrix::allocateData(bool Complex){
    freeData();
    conjugated=false;
    if(shape.total==0)return;
    // NOTE: all matrices have 0 as the last value at data+total
    if(Complex){cdata=new complex<double>[shape.total];cornerData(cdata);*(cdata+shape.total-1)=0.;}
    else       {rdata=new          double[shape.total];cornerData(rdata);*(rdata+shape.total-1)=0.;}
}

UseMatrix::UseMatrix(const UseMatrix & other) :
    rdata(0),cdata(0),conjugated(false),symm(new Symmetry(other.symm->symm))
{
    trf=new TriFactorLU();
    shape=Shape(other.rows(),other.cols(),other.shape.lowBw-1,other.shape.upBw-1);
    allocateData(other.cdata!=0);

    if(cdata!=0)ternaryData(BinaryInAssign,1.,cdata,other,1.,other.cdata,*this,1.,cdata);
    else        ternaryData(BinaryInAssign,1.,rdata,other,1.,other.rdata,*this,1.,rdata);
}


void UseMatrix::swap(UseMatrix & B){
    std::swap(shape,B.shape); // note: we would need a specialized swap here
    std::swap(rdata,B.rdata);
    std::swap(cdata,B.cdata);
    std::swap(conjugated,B.conjugated);
    std::swap(trf,B.trf);
}

UseMatrix & UseMatrix::assignMatrix(const UseMatrix & other, symmetry Symm, bool ForceFull){
    if(Symm and !other.isSymmetric())ABORT("matrix is not symmetric, cannot cast into symmetric storage");
    if(&other!=this){
        delete trf; trf=0;
        conjugated=other.conjugated;
        conjugated=false;
        if(isBasic() and not sharesData(other)){
            // re-allocate in the correct shape
            if(ForceFull)shape=Shape(other.shape.nrows,other.shape.ncols);
            else shape=Shape(other.shape.nrows,other.shape.ncols,other.shape.lowBw-1,other.shape.upBw-1,Symm);
            allocateData(other.cdata!=0);
        } else {
            //            ABORT("case not handled yet");
            if(rows()!=other.rows() or cols()!=other.cols())
                ABORT("matrix dimensions do not match - cannot assign\n   "+shape.str()+"\n!="+other.shape.str());
        }

        if     (cdata!=0)ternaryData(BinaryInAssign,1.,cdata,other,1.,other.cdata,*this,1.,cdata);
        else if(rdata!=0)ternaryData(BinaryInAssign,1.,rdata,other,1.,other.rdata,*this,1.,rdata);
    }
    return *this;
}

UseMatrix::EigenMethod UseMatrix::eigenMethod=UseMatrix::automatic;
void UseMatrix::eigenValues(UseMatrix &Val, const UseMatrix &Ovr) const{
    UseMatrix dum(1,1);
    eigen(Val,dum,Ovr,false);
}

void UseMatrix::eigenBlock(UseMatrix &Val, UseMatrix &Vec, UseMatrix &Ovr, bool Vectors)
{
    UseMatrix dum;
    eigenBlock(Val,Ovr,Vectors,Vec,false,dum);
}

void UseMatrix::eigenBlock(UseMatrix &Val, UseMatrix &Ovr, bool RightVectors, UseMatrix &RVec, bool DualVectors, UseMatrix &DVec){

    // get blocking
    vector<vector<unsigned int> > blocks=blocking(vector<const UseMatrix*>(1,&Ovr),1.e-12);
    PrintOutput::DEVmessage("eigenBlock found independent subblocks: "+tools::str(int(blocks.size())));

    if(rows()<smallDimensions){
        // do not block-decompose small matrices
        eigen(Val,Ovr,RightVectors,RVec,DualVectors,DVec);
    }

    else {

        unsigned int iVals=0; // block starting index
        // loop through blocks
        Val=UseMatrix::Zero(cols(),1);
        RVec=UseMatrix();
        DVec=UseMatrix();
        if(RightVectors)RVec=UseMatrix::Zero(cols(),cols());
        if( DualVectors)DVec=UseMatrix::Zero(cols(),cols());
        for(unsigned int k=0;k<blocks.size();k++){

            // extract matrix and overlap
            UseMatrix bMat=UseMatrix::Zero(blocks[k].size(),blocks[k].size());
            UseMatrix bOvr=UseMatrix::Zero(blocks[k].size(),blocks[k].size());
            UseMatrix bVal(bMat.cols(),1),bRVec,bLVec;
            for (unsigned int i=0;i<blocks[k].size();i++){
                for (unsigned int j=0;j<blocks[k].size();j++){
                    bMat(i,j)=operator()(blocks[k][i],blocks[k][j]);
                    bOvr(i,j)=Ovr(blocks[k][i],blocks[k][j]);
                }
            }
            // solve sub-block eigenproblem
            bMat.eigen(bVal,bOvr,RightVectors,bRVec,DualVectors,bLVec);
            if(RightVectors and not DualVectors){
                if(not eigenOrthonormalize(bVal,bRVec,bOvr))
                    PrintOutput::warning("not a symmetric or hermitian problem, cannot (pseudo-)orthonormalize",5);
            }

            // place block eigenvalues and vectors into global positions
            Val.block(iVals,0,bVal.rows(),1)=bVal;
            for(unsigned int i=0;i<bRVec.cols();i++){
                if(RightVectors)RVec.block(blocks[k][i],iVals,1,bRVec.cols())=bRVec.row(i);
                if( DualVectors)DVec.block(blocks[k][i],iVals,1,bRVec.cols())=bLVec.row(i);
            }
            iVals+=bVal.size();
        }
    }
}

/// return true if sucessfully orthonormalized or pseudo-orthonormalized
bool UseMatrix::eigenOrthonormalize(UseMatrix &Val, UseMatrix &Vec, const UseMatrix &Ovr, double Eps) const {
    // find blocks of near-degenerated values
    vector<bool> use(Val.size(),true);

    bool pseudo=not(isHermitian() and Ovr.isHermitian());
    if(pseudo and not(isSymmetric() and Ovr.isSymmetric()))return false;

    vector<unsigned int>block;
    do {
        block.clear();
        for(unsigned int k=0;k<Val.size();k++){
            if(use[k] and (block.size()==0 or
                           abs(Val(block[0])-Val(k))<Eps*max(1.,abs(Val(k))))){
                block.push_back(k);
                use[k]=false;
            }
        }
        // Schmidt-orthonormalize within block (may be single function!)
        if(block.size()>0)Vec.gramSchmidt(Ovr,block,pseudo);

        // advance to next block
    } while(block.size()>0);
    return true;
}

/// considering the matrix as a block of eigenvectors, (pseudo)-orthonormalizes Subset wrt Ovr
void UseMatrix::gramSchmidt(const UseMatrix & Ovr,vector<unsigned int> Subset,bool Pseudo){
    if(Subset.size()==0)
        for(unsigned int k=0;k<cols();k++)Subset.push_back(k);
    UseMatrix OvrI(rows(),1);
    for (unsigned int i=0;i<Subset.size();i++){
        if(Ovr.size()!=0)OvrI=Ovr*col(Subset[i]);
        else             OvrI=col(Subset[i]);
        complex<double>ovrij;
        for(unsigned int j=0;j<i;j++){
            if(Pseudo) ovrij=(col(Subset[j]).transpose()*OvrI)(0,0);
            else       ovrij=(col(Subset[j]).adjoint()*OvrI)(0,0);
            col(Subset[i])-=col(Subset[j])*ovrij;
        }
        if(Pseudo) ovrij=(col(Subset[i]).transpose()*OvrI)(0,0);
        else       ovrij=(col(Subset[i]).adjoint()*OvrI)(0,0);
        if(abs(ovrij)<1.e-28)ABORT(Str("vectors are (pseudo-)linearly dependent")+ovrij+"pseudo="+Pseudo);
        col(Subset[i])*=1./sqrt(ovrij);
    }
}
void UseMatrix::eigen(UseMatrix &Val, UseMatrix &Vec, const UseMatrix &Ovr, bool Vectors) const
{
    UseMatrix dum;
    eigen(Val,Ovr,Vectors,Vec,false,dum);
}

std::string sortDuals(Eigen::MatrixXcd Dval,Eigen::MatrixXcd& Dvec,const Eigen::MatrixXcd Rval,const Eigen::MatrixXcd Rvec){
    double Eps=1e-12;
    std::vector<std::vector<int>> rDeg;

    // make sure eigenvalues match
    for(int k=0;k<Rval.size();k++){
        // collect the degenerate eigenvalues
        bool found=false;
        for(auto &n: rDeg)
            if(std::abs(Rval(n[0])-Rval(k))<Eps){
                n.push_back(k);
                found=true;
            }
        if(not found)rDeg.push_back(std::vector<int>(1,k));
    }

    std::vector<std::vector<int>> dDeg(rDeg.size());
    for(int l=0;l<Dval.size();l++){
        bool found=false;
        for(size_t n=0;n<rDeg.size();n++){
            if(std::abs(Rval(rDeg[n][0])-std::conj(Dval(l)))<Eps*std::max(1.,std::abs(Dval(l)))){
                dDeg[n].push_back(l);
                found=true;
            }
        }
        if(not found)DEVABORT(Sstr+"lhs eigenvalue not in rhs eigenvalues:\n"
                              +Dval.transpose()+"\n"+Rval.transpose());
    }

    // sort for Duals for maximal overlap within degenerate
    for(size_t n=0;n<rDeg.size();n++){
        for(auto k: rDeg[n]){
            double maxOvr=0;
            int match=-1;
            for(int l=0;l<Dval.size();l++){
                if(std::abs(Rval(k)-std::conj(Dval(k)))<Eps){
                    std::complex<double>ovr(0.);
                    for(int r=0;r<Dvec.rows();r++)ovr+=std::conj(Dvec(r,l))*Rvec(r,k);
                    if(std::abs(ovr)>maxOvr){
                        maxOvr=std::abs(ovr);
                        match=l;
                    }
                }
            }
            if(match==-1){
                DEVABORT(Sstr+"left and right eigenvalues do not seem to match"+dDeg[n]
                         +"Vec:"+k+Rvec.col(k).transpose()
                         +"\n (unclarified problem with lhs eigenvectors");
            }
            Dvec.col(k).swap(Dvec.col(match));
            std::swap(Dval(k),Dval(match));
        }
    }
    DEVABORT("return value missing");
}


void UseMatrix::eigen(UseMatrix &Val, const UseMatrix &Ovr, bool RightVectors, UseMatrix &RVec, bool DualVectors, UseMatrix &DVec) const {

    // special cases
    if(cols()*rows()==0){
        ABORT("zero size matrix - eigenvectors not defined");
    } else if (cols()*rows()==1){
        complex<double> ovr=Ovr(0,0).complex();
        if(ovr==0.)ABORT("1x1 overlap = 0");
        Val=*this;
        Val/=ovr;
        if(RightVectors)RVec=UseMatrix::Constant(1,1,1./sqrt(ovr));
        if(DualVectors)DVec=RVec;
        return;
    }

    STARTDEBUG(eigen);
    // decision tree: try exploiting symmetry and band structure
    EigenMethod Method=UseMatrix::eigenMethod;
    if(Method==automatic){
        Method=general;
        if(not RightVectors and isSymmetric()){ // note: the symmetric banded solver does not do vectors at present
            if(Ovr.isSymmetric()){
                int trueBand=1+max(trueSubD(),Ovr.trueSubD())+max(trueSuperD(),Ovr.trueSuperD());
                if(trueBand*trueBand*cols()*2<cols()*cols()*cols())Method=symmetric_banded;
                else Method=symmetric_full;
            }    }
    }
    Val=UseMatrix::Zero(rows(),1);

    switch (Method){
    case symmetric_full:
    case general:
    {
        UseMatrix mat(*this);
        if(Ovr.size()!=0){
            Ovr.solve(mat);
        }

        char jobv='n',jobl='n';
        if(RightVectors){
            RVec=UseMatrix::Zero(rows(),cols());
            jobv='v';
        }
        if(DualVectors){
            DVec=UseMatrix::Zero(rows(),cols());
            jobl='v';
        }
        Val=UseMatrix::Zero(rows(),1);

        int info,ilo,ihi;
        double abnrm;
        vector<double> scale(cols()),rconde(cols()),rcondv(cols());
        const char balance='N'; //CAUTION: attemps to improve by choosing ='B' proved counter-productive
//#ifndef _NO_LAPACKE_
//        info=LAPACKE_zgeevx(LAPACK_COL_MAJOR,balance,jobl,jobv,'N',cols(),
//                            mat.data(0,0),mat.leadDim(),Val.data(),DVec.data(),DVec.leadDim(),RVec.data(0,0),RVec.leadDim(),
//                            &ilo,&ihi,scale.data(),&abnrm,rconde.data(),rcondv.data());

//        if(info!=0)ABORT("zgeev failed, info="+tools::str(info));
//        for(size_t k=0;k<rcondv.size();k++){
//            if(rconde[k]>1.e8)PrintOutput::DEVwarning(Str("large eigenvalue  condition number")+rconde[k]);
//            if(rcondv[k]>1.e8)PrintOutput::DEVwarning(Str("large eigenvector condition number")+rcondv[k]);
//        }
//#else

        if(mat.leadDim()!=mat.rows())DEVABORT("must have leadDim==rows");
        Eigen::VectorXcd eval;
        if(mat.isHermitian()){
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>slv(Eigen::Map<Eigen::MatrixXcd>(mat.data(),mat.rows(),mat.cols()),
                                                               RightVectors?Eigen::ComputeEigenvectors:Eigen::EigenvaluesOnly);
            Eigen::Map<Eigen::VectorXcd>(Val.data(),Val.size())=slv.eigenvalues();
            if(RightVectors)Eigen::Map<Eigen::MatrixXcd>(RVec.data(),RVec.rows(),RVec.cols())=slv.eigenvectors();
            eval=slv.eigenvalues();
        }
        else {
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> slv(Eigen::Map<Eigen::MatrixXcd>(mat.data(),int(mat.rows()),int(mat.cols())),RightVectors);
            eval=slv.eigenvalues();
            Eigen::Map<Eigen::VectorXcd>(Val.data(),Val.size())=slv.eigenvalues();
            if(RightVectors)Eigen::Map<Eigen::MatrixXcd>(RVec.data(),RVec.rows(),RVec.cols())=slv.eigenvectors();
        }

        if(DualVectors){
            if(mat.isHermitian())
                DVec=RVec;
            else if(mat.isSymmetric(1e-12))
                for(size_t k=0;k<DVec.size();k++)DVec.data()[k]=conj(RVec.data()[k]);
            else {
                Eigen::Map<Eigen::MatrixXcd>(mat.data(),int(mat.rows()),int(mat.cols())).adjointInPlace();
                Eigen::ComplexEigenSolver<Eigen::MatrixXcd> slv(Eigen::Map<Eigen::MatrixXcd>(mat.data(),int(mat.rows()),int(mat.cols())),true);
                for(size_t k=0;k<Val.size();k++)
                    if(std::abs(Val(k)-std::conj(slv.eigenvalues()(k)))>1.e-12*std::max(1.,std::abs(Val(k).complex()))){
                        Sstr+"inaccurate eigenvalues"+k+Val(k)+std::conj(slv.eigenvalues()(k))+Sendl;
                    }

                Eigen::MatrixXcd dvec=slv.eigenvectors();
                sortDuals(slv.eigenvalues(),dvec,eval,Eigen::Map<Eigen::MatrixXcd>(RVec.data(),RVec.rows(),RVec.cols()));
                Eigen::Map<Eigen::MatrixXcd>(DVec.data(),DVec.rows(),DVec.cols())=dvec;
            }
        }
//#endif
        if(DualVectors)
            for(size_t k=0;k<DVec.size();k++)DVec.data()[k]=conj(DVec.data()[k]);
        if(DualVectors){
            for(size_t k=0;k<RVec.cols();k++){
                double ovrMax=0;
                UseMatrix rCol=RVec.block(0,k,RVec.rows(),1);
                for(size_t l=0;l<DVec.cols();l++){
                    std::complex<double> ovr=0;
                    for(size_t r=0;r<DVec.rows();r++)
                        ovr+=DVec.data()[r+l*DVec.rows()]*RVec.data()[r+k*RVec.rows()];
                    if(std::abs(ovr)>ovrMax){
                        ovrMax=std::abs(ovr);
                    }
                }
            }
        }
    }
        break;
    case symmetric_banded:
    {
        Sstr+"symmetric banded solver"+Sendl;
        // re-store matrices in symmetric banded format
        if(Ovr.size()!=size())ABORT("for banded, Ovr default not implemented");
        UseMatrix ab,bb;
        STARTDEBUG(compress);
        ab.assignMatrix(band(trueSubD(),trueSuperD()),symmetric);
        bb.assignMatrix(Ovr.band(Ovr.trueSubD(),Ovr.trueSuperD()),symmetric);
        STOPDEBUG(compress);
        ab.printShape("banded");
        bb.printShape("overlap");


        if(RightVectors or DualVectors)ABORT("banded symmetric solver cannot do eigenvectors at present");
#ifdef _USE_ZSLAPACK_
        zslapack::bgv(RightVectors,ab.cols(),ab.subD(),bb.subD(),ab.cdata,ab.leadDim(),bb.cdata,bb.leadDim(),Val.cdata,0,1);
#else
        ABORT("for banded complex symmetric, compile with -D_USE_ZSLAPACK_");
#endif
    }
        break;

    default:
        ABORT("eigensolver not implemented for this case");
    }
    STOPDEBUG(eigen);

    // reset eigenMethod to automatic
    UseMatrix::eigenMethod=UseMatrix::automatic;
}

//NOTE: this algorithm has been templated into tRecX::matrixBlocking (remove from here eventually)
vector<vector<unsigned int> > UseMatrix::blocking(const std::vector<const UseMatrix *> &Mat, double Eps) const{
    // check dimensions
    if(rows()!=cols())ABORT("implemented only for square matrices");
    for(unsigned int k=0;k<Mat.size();k++){
        if(rows()!=Mat[k]->rows())ABORT(Str("rows do not match")+rows()+Mat[k]->rows());
        if(cols()!=Mat[k]->cols())ABORT(Str("columns do not match")+cols()+Mat[k]->cols());
    }

    // list of block
    vector<vector<unsigned int> > blocks(0);

    // find connected subsets
    unsigned int next=0;
    vector<bool> Use(rows(),true);
    while (next<rows()) {
        blocks.push_back(vector<unsigned int>(1,next));
        Use[next]=false;
        connect(next,blocks.back(),Use,Mat,Eps);

        // find next unused
        for(;next<Use.size();next++)
            if(Use[next])
                break;
    }
    //    for(unsigned int k=0;k<blocks.size();k++)cout<<"block: "<<tools::str(blocks[k])<<endl;
    return blocks;
}

void UseMatrix::connect(unsigned int From, vector<unsigned int> & Conns, vector<bool> & Use,
                        std::vector<const UseMatrix *> Other, double Eps) const {
    if(Use.size()==0)Use.assign(rows(),true);
    if(Use.size()!=rows())ABORT("Use.size() must match rows()");
    Other.push_back(this);

    // get all indices that connect From and that are not taken yet
    vector<unsigned int> newConns;
    for(unsigned int k=0;k<rows();k++){
        if(Use[k]){

            // check whether k connects to present index
            bool connected=k==From;
            for(unsigned int i=0;i<Other.size();i++){
                if(connected)break;
                connected=connected or
                        pow(abs(Other[i]->operator()(k,From)),2)
                        >Eps*Eps*abs(Other[i]->operator()(k,k)*Other[i]->operator()(From,From));
            }

            if(connected){
                Use[k]=false;
                Conns.push_back(k);
                if(k!=From)newConns.push_back(k);
            }
        }
    }
    Other.pop_back();

    // add all connections from all new indices
    for (unsigned int k=0;k<newConns.size();k++)
        connect(newConns[k],Conns,Use,Other,Eps);
}



// construct column vector from data array
UseMatrix UseMatrix::UseMap(complex<double>* Data, int Nrows) {return UseMap(Data,Nrows,1);}

// construct full matrix from data array
UseMatrix UseMatrix::UseMap(std::complex<double> *Data, int Nrows, int Ncols) {
    UseMatrix C;
    C.shape=Shape(Nrows,Ncols);
    C.shape.locationIJ=full_normal_map; // will be treated as derived matrix
    C.conjugated=false;
    C.rdata=0;
    if(Nrows==0 or Ncols==0)C.cdata=0;
    else                    C.cdata=Data;
    //    C.isIdentityConst=C.isIdentity();
    return C;
}

// default constructor
UseMatrix::UseMatrix()
    :rdata(0),cdata(0),conjugated(false),symm(new Symmetry(unknown))
{trf=new TriFactorLU();shape=Shape();allocateData(true);}

UseMatrix::UseMatrix(unsigned int Nrows,unsigned int Ncols) :
    rdata(0),cdata(0),conjugated(false),symm(new Symmetry(unknown))
{
    trf=new TriFactorLU();
    shape=Shape(Nrows,Ncols);
    allocateData(true);
    //    isIdentityConst = isIdentity();
}
UseMatrix::UseMatrix(unsigned int Nrows,unsigned int Ncols,unsigned int SubDiag,unsigned int SuperDiag, symmetry Symm) :
    rdata(0),cdata(0),conjugated(false),symm(new Symmetry(unknown))
{
    trf=new TriFactorLU();
    if(SubDiag==Nrows-1 or SuperDiag==Ncols-1)shape=Shape(Nrows,Ncols); // is actually a full matrix
    else                                      shape=Shape(Nrows,Ncols,SubDiag,SuperDiag,Symm);
    allocateData(true);
    //    isIdentityConst = isIdentity();
}

// return Data containing pointers to storage location
UseMatrix::Data UseMatrix::operator()(unsigned int I) const {
    if(shape.ncols==1)return operator()(I,0);
    else if(shape.nrows==1)return operator()(0,I);
    else ABORT("single index only for single row/single column matrices");
    abort();
}
UseMatrix::Data UseMatrix::operator()(unsigned int I, unsigned int J) const {
    if(I<0        or J<0       )ABORT("matrix subscripts < 0");
    if(size()==0)ABORT("zero size matrix - cannot access elements");
    if(I>rows()-1 or J>cols()-1){
        cout<<"I,J "<<I<<", "<<J<<endl;
        printShape("");
        ABORT("matrix subscripts >= rows/cols");
    }
    return Data(this,I,J);
}

UseMatrix::Data::Data(const UseMatrix * M, int I, int J):m(M) {
    unsigned int loc=M->shape.locIJ(I,J);
    conjugated=M->conjugated;
    if(M->cdata!=0){
        cZero=M->cdata+M->shape.total-1;
        c=M->cdata+loc;
        r=0;
    } else {
        rZero=M->rdata+M->shape.total-1;
        c=0;
        r=M->rdata+loc;
    }
}

// check whether data can be written to and reset factorization
void UseMatrix::Data::changeData(){
    delete(m->trf);const_cast<UseMatrix*>(m)->trf=0;
    if((c!=0 and c==cZero) or (r!=0 and r==rZero))
        ABORT("cannot write to matrix - element outside range (band, block)");
}

void UseMatrix::Data::operator=(double D){
    changeData();
    if     (c!=0){(*c)=D;}
    else if(r!=0){(*r)=D;}
    else ABORT("UseMatrix::Data empty, cannot assign value");
}

void UseMatrix::Data::operator=(std::complex<double> C){
    changeData();
    if(c!=0){
        if(conjugated)*c=conj(C);
        else          *c=C;
    }
    else if(r!=0 and std::imag(C)==0.)*r=std::real(C);
    else ABORT("cannot insert complex value into real matrix");
}
void UseMatrix::Data::operator*=(std::complex<double> C){
    changeData();
    if(c!=0){
        if(conjugated)*c*=conj(C);
        else          *c*=C;
    }
    else if(r!=0 and std::imag(C)==0.)*r*=std::real(C);
    else ABORT("cannot insert complex value into real matrix");
}
void UseMatrix::Data::operator/=(std::complex<double> C){
    changeData();
    if(c!=0){
        if(conjugated)*c/=conj(C);
        else          *c/=C;
    }
    else if(r!=0 and std::imag(C)==0.)*r/=std::real(C);
    else ABORT("cannot insert complex value into real matrix");
}
void UseMatrix::Data::operator+=(std::complex<double> C){
    changeData();
    if(c!=0){
        if(conjugated)*c+=conj(C);
        else          *c+=C;
    }
    else if(r!=0 and std::imag(C)==0.)*r+=std::real(C);
    else ABORT("cannot insert complex value into real matrix");
}
void UseMatrix::Data::operator-=(std::complex<double> C){
    changeData();
    if(c!=0){
        if(conjugated)*c-=conj(C);
        else          *c-=C;
    }
    else if(r!=0 and std::imag(C)==0.)*r-=std::real(C);
    else ABORT("cannot insert complex value into real matrix");
}

void UseMatrix::Data::operator=(Data Dat){
    changeData();
    if(c!=0){
        if(Dat.c!=0){
            if(conjugated!=Dat.conjugated)*c=conj(*Dat.c);
            else                          *c=*Dat.c;
        }
        else if(Dat.r!=0)*r=*(Dat.r);
        else ABORT("UseMatrix::Data emtpy, cannot extract value");
    } else if(r!=0){
        if     (Dat.r!=0)*r=*(Dat.r);
        else if(Dat.c!=0 and std::imag(*(Dat.c))==0.)*r=std::real(*(Dat.c));
        ABORT("cannot insert complex value into real matrix");
    } else {
        if(Dat.c!=0 or Dat.r!=0)ABORT("UseMatrix:Data empty, cannot assign to it");
    }
}

bool UseMatrix::Data::operator==(const Data & Dat) const {
    if(c!=0){
        if(Dat.c!=0){
            if(conjugated!=Dat.conjugated)return *c==conj(*Dat.c);
            else                          return *c==*Dat.c;
        }
        else if(Dat.r!=0 and std::imag(*c)==0.)return std::real(*c)==*(Dat.r);
        else return false;
    } else if (r!=0){
        if     (Dat.r!=0)return *r==*(Dat.r);
        else if(Dat.c!=0 and std::imag(*(Dat.c))==0.)return std::real(*(Dat.c))==*r;
        else return false;
    } else {
        return (Dat.c==0) and (Dat.r==0); // both emtpy
    }
}

#ifdef _CBLAS_
CBLAS_TRANSPOSE UseMatrix::cblasTrans() const {
    if(shape.storeTranspose()){
        if(conjugated)return CblasConjTrans;
        return CblasTrans;
        //        if(conjugated)return CblasConjTrans;
        //        return CblasTrans;
    }
    else return CblasNoTrans;
}
#endif

bool UseMatrix::isFull() const {
    if(isShape(full_normal))return true;
    if(shape.stored!=0)return shape.stored->storeFull();
    return false;
}

bool UseMatrix::isBand() const {return shape.storeBand();}
bool UseMatrix::isBandNormal() const {return (shape.locationIJ==band_normal);}
bool UseMatrix::Shape::storeBand() const {
    if(stored!=0)return stored->storeBand();
    if(locationIJ==band_normal)return true;
    if(locationIJ==band_symmetric)return true;
    return false;
}

UseMatrix & UseMatrix::operator*=(const UseMatrix & B){
    UseMatrix C(rows(),B.cols());
    multiply(B,C,1.,0.);
    swap(C);
    return *this;
}


// C=c*C+a*A*B
void UseMatrix::multiply(const UseMatrix & B, UseMatrix & C, const z & a, const z & c) const
{
    if(rows()==0 or cols()==0 or B.cols()==0){
        C=UseMatrix();
        return;
    }

    if(cols()!=B.rows()){
        shape.print(  "A*B: A");
        B.shape.print("     B");
        ABORT("cannot multiply matrices - dimensions do not match");
    }

    // if matrices share data, need a copy
    const UseMatrix *pA=this,*pB=&B;
    if(sharesData(C)){pA=new UseMatrix(*this);trf->reset();} // C=this will be modified, reset factorization
    if(sharesData(B))pB=new UseMatrix(B);

    if(rdata!=0)ABORT("real not implemented yet");

    // array start
    if(isFull() and B.isFull()){
        unsigned int na=pA->shape.locIJ(0,0);
        unsigned int nb=pB->shape.locIJ(0,0);
        unsigned int nc=  C.shape.locIJ(0,0);

        // do the BLAS operation
#ifdef _CBLAS_
        cblas_zgemm(CblasColMajor,pA->cblasTrans(),pB->cblasTrans(),
                    pA->rows(),pB->cols(),pA->cols(),
                    reinterpret_cast<const double*>(&a),reinterpret_cast<const double*>(pA->cdata+na),
                    pA->shape.leadDim,reinterpret_cast<double*>(pB->cdata+nb), pB->shape.leadDim,
                    reinterpret_cast<const double*>(&c), reinterpret_cast<double*>(C.cdata+nc),
                    C.shape.leadDim);
#else
        Eigen::Map<Eigen::MatrixXcd>(C.cdata+nc,C.leadDim(),C.cols())*=c;
        if(not pA->shape.storeTranspose() and not pA->conjugated and not pA->shape.storeTranspose() and not pA->conjugated){
            Eigen::Map<Eigen::MatrixXcd>(C.cdata+nc,C.leadDim(),C.cols())+=a
                    *Eigen::Map<Eigen::MatrixXcd>(pA->cdata+na,pA->leadDim(),pA->cols())
                    *Eigen::Map<Eigen::MatrixXcd>(pB->cdata+nb,pB->leadDim(),pB->cols());
        }
        else if(pA->shape.storeTranspose() and not pA->conjugated and  pA->shape.storeTranspose() and not pA->conjugated){
            Eigen::Map<Eigen::MatrixXcd>(C.cdata+nc,C.leadDim(),C.cols())+=a
                    *Eigen::Map<Eigen::MatrixXcd>(pA->cdata+na,pA->leadDim(),pA->rows()).transpose()
                    *Eigen::Map<Eigen::MatrixXcd>(pB->cdata+nb,pB->leadDim(),pB->rows()).transpose();
        }
        else
            DEVABORT(Sstr+"no CBLAS available - cannot use UseMatrix::multiply for full"
                     +pA->shape.storeTranspose()+pA->conjugated+pA->shape.storeTranspose()+pA->conjugated);
#endif
    }
    else if(isShape(band_normal) and B.isFull()){


        // series of matrix-vector products
#ifdef _CBLAS_
        int incB=1,incC=1;
        if(C.shape.storeTranspose())incC=C.leadDim();
        if(B.shape.storeTranspose())incB=B.leadDim();
        complex<double>*A0=pA->cdata+pA->shape.row0+pA->shape.col0*pA->shape.leadDim;
        for(unsigned int col=0;col<B.cols();col++)
            cblas_zgbmv(CblasColMajor,pA->cblasTrans(),rows(),cols(),subD(),superD(),
                        reinterpret_cast<const double*>(&a),reinterpret_cast<const double*>(A0),pA->shape.leadDim,
                        reinterpret_cast<double*>(pB->cdata+pB->shape.locIJ(0,col)),incB,
                        reinterpret_cast<const double*>(&c),reinterpret_cast<double*>(C.cdata+C.shape.locIJ(0,col)),
                        incC);
#else
        DEVABORT("no CBLAS available - cannot use UseMatrix::multiply for banded");
#endif
    }
    else {
        ABORT("not implemented for shapes "+shape.strLocIJ()+"*"+B.shape.strLocIJ());
    }

    if(pA!=this)delete pA;
    if(pB!=&B  )delete pB;
}

// all go through ternary ... can be advantegously split....
void UseMatrix::scalarOp(Operation Op, std::complex<double> a){ternaryOp(Op,*this,*this,true,1.,1.,a);}
void UseMatrix::binaryOp(Operation Op, const UseMatrix & B, std::complex<double> a,std::complex<double> b){ternaryOp(Op,B,*this,true);}

// ternary matrix operations for all shapes
void UseMatrix::ternaryOp(Operation Op, const UseMatrix & B, UseMatrix & C, bool KeepC, const z a,const z b,const z c) const {

    // the matrix will be modified, its factorization becomes invalid
    if(trf!=0 and Op.inPlace)trf->reset();

    if(Op.term==ternary and (rows()!=B.rows() or cols()!=B.cols())){
        shape.print("A");
        B.shape.print("B");
        C.shape.print("C");
        ABORT("cannot combine - dimensions do not match");
    }

    // zero size matrix - no action taken;
    if(size()==0)return;

    // determine band-width of output storage
    unsigned int subOut=0,superOut=0;
    switch (Op.kind) {
    case assign:
    case add:
    case subtract:
        subOut  =max(shape.lowBw,B.shape.lowBw)-1;
        superOut=max(shape.upBw, B.shape.upBw)-1;
        break;
    case product:
        subOut  =min(shape.lowBw,B.shape.lowBw)-1;
        superOut=min(shape.upBw, B.shape.upBw)-1;
        break;
    case quotient:
        subOut  =shape.lowBw-1;
        superOut=shape.upBw-1;
        if(subOut+1<B.shape.lowBw or superOut+1<B.shape.upBw)ABORT("for quotient, second matrix must be at least as wide as first");
        break;
    default:
        ABORT("undefined ternary operator");
    }

    // re-shape input if needed
    UseMatrix rA,rB;
    rA=  reband(subOut,superOut);
    rB=B.reband(subOut,superOut);

    // prepare output
    if(Op.inPlace or KeepC){
        // keep data of C (either by parameter: KeepC or by function type: Op.inPlace)
        if(rows()!=C.rows() or cols()!=C.cols())ABORT("cannot keep output - output dimensions do not match");
        C=C.reband(subOut,superOut);
    } else {
        // supersede with new basic matrix
        C=UseMatrix(rows(),B.cols(),subOut,superOut);
    }

    // dynamic type selection
    const int Za=1,Zb=2,Zc=4,za=8,zb=16,zc=32;
    const int RRRrrr=0,RRZrrr=Zc,ZRZrrr=Za+Zc,RZZrrr=Zb+Zc,ZZZrrr=Za+Zb+Zc;
    const int RRZzrr=RRZrrr+za; // add more cases as needed
    const int ZRZzrr=ZRZrrr+za; // add more cases as needed
    const int ZRZrrz=ZRZrrr+zc; // add more cases as needed
    const int ZRZrzz=ZRZrrz+zb; // add more cases as needed
    const int ZZZrrz=ZZZrrr+zc; // add more cases as needed

    int code=0;
    if(  cdata!=0 )code+=1;
    if(B.cdata!=0 )code+=2;
    if(C.cdata!=0 )code+=4;
    if(imag(a)!=0.)code+=8;
    if(imag(b)!=0.)code+=16;
    if(imag(c)!=0.)code+=32;

    // apply the operations
    switch (code) {
    case ZZZrrr: rA.ternaryData(Op,real(a),rA.cdata,rB,real(b),rB.cdata,C,real(c),C.cdata);break;
    case ZZZrrz: rA.ternaryData(Op,real(a),rA.cdata,rB,real(b),rB.cdata,C,     c ,C.cdata);break;

    case RZZrrr: rA.ternaryData(Op,real(a),rA.rdata,rB,real(b),rB.cdata,C,real(c),C.cdata);break;

    case ZRZrrr: rA.ternaryData(Op,real(a),rA.cdata,rB,real(b),rB.rdata,C,real(c),C.cdata);break;
    case ZRZzrr: rA.ternaryData(Op,     a ,rA.cdata,rB,real(b),rB.rdata,C,real(c),C.cdata);break;
    case ZRZrzz: rA.ternaryData(Op,real(a),rA.cdata,rB,     b ,rB.rdata,C,     c ,C.cdata);break;

    case RRZrrr: rA.ternaryData(Op,real(a),rA.rdata,rB,real(b),rB.rdata,C,real(c),C.cdata);break;
    case RRZzrr: rA.ternaryData(Op,     a ,rA.rdata,rB,real(b),rB.rdata,C,real(c),C.cdata);break;

    case RRRrrr: rA.ternaryData(Op,real(a),rA.rdata,rB,real(b),rB.rdata,C,real(c),C.rdata);break;

    default: ABORT("combination of data types not implemented: ");
    }
}

template<class Ad,class Bd,class Cd,class As,class Bs,class Cs>
void UseMatrix::ternaryData(Operation Op,        const As as, const Ad & Adata,
                            const UseMatrix & B, const Bs bs, const Bd & Bdata,
                            const UseMatrix & C, const Cs cs, const Cd & Cdata
                            ) const {
    // ternary element-wise operations on matrix data for all shapes
    //
    // C<-f(A,B,C)
    //
    // core routine of the matrix class

    if((shape==B.shape) and (shape==C.shape) and shape.storeContiguous()){
        // contiguous data on identical shapes

        // for divisions, make sure there are no zeros in the unused corners of the storage
        if(Op.kind==quotient)B.cornerData(Bdata);

        // beginning and end of iterators
        Ad a=Adata+shape.locIJ(0,0),Aend=Adata+  shape.total;
        Bd b=Bdata+shape.locIJ(0,0),Bend=Bdata+B.shape.total;
        Cd c=Cdata+shape.locIJ(0,0),Cend=Cdata+C.shape.total;

        rangeOp(Op.kind,Op.term,Op.inPlace,C.conjugated!=conjugated,C.conjugated!=B.conjugated,
                a,Aend,1,b,Bend,1,c,Cend,1,as,bs,cs);
    } else {
        // column-wise operations

        // get increments (depending on transposition)


        int Ainc=storeIncI();
        int Binc=B.storeIncI();
        int Cinc=C.storeIncI();

        // check band-width
        if(  (Op.term==ternary and (  shape.lowBw>C.shape.lowBw or   shape.upBw>C.shape.upBw)) or
             (Op.term!=scalar  and (B.shape.lowBw>C.shape.lowBw or B.shape.upBw>C.shape.upBw)))
        {
            shape.print(  "A");
            B.shape.print("B");
            C.shape.print("C");
            ABORT("In f(A,B)->C: band width too small");
        }

        // for supplementing when opperating with symmetric storage
        bool hasSymm=shape.storeSymmetric() or B.shape.storeSymmetric() or C.shape.storeSymmetric();
        int AincT=1,BincT=1,CincT=1;
        if(hasSymm and not C.shape.storeSymmetric()){
            if(!  shape.storeTranspose() and !shape.storeSymmetric())AincT=  storeIncI();
            if(!B.shape.storeTranspose() and !shape.storeSymmetric())BincT=B.storeIncI();
            if(!C.shape.storeTranspose() and !shape.storeSymmetric())CincT=C.storeIncI();
        }

        for (unsigned int j=0;j<C.cols();j++){

            // get starting locations in the column
            unsigned int imin=max(     0,(int) j- int(C.shape.upBw)+1);
            if(hasSymm)imin=j;

            Ad a= Adata+  shape.locIJ(imin,j);
            Bd b= Bdata+B.shape.locIJ(imin,j);
            Cd c= Cdata+C.shape.locIJ(imin,j);

            // get end for iterator
            unsigned int imax=min((int) rows()-1,(int) j+ int(C.shape.lowBw)-1);
            Ad Aend=Adata+  shape.locIJ(imax,j)+1;
            Bd Bend=Bdata+B.shape.locIJ(imax,j)+1;
            Cd Cend=Cdata+C.shape.locIJ(imax,j)+1;

            rangeOp(Op.kind,Op.term,Op.inPlace,C.conjugated!=conjugated,C.conjugated!=B.conjugated,
                    a,Aend,Ainc,b,Bend,Binc,c,Cend,Cinc,as,bs,cs);

            if(hasSymm and not C.shape.storeSymmetric()){
                // operate on symmetric part if needed

                imin=max(     0,(int) j- int(C.shape.upBw)+1);
                a   =Adata+  shape.locIJ(imin,j);
                b   =Bdata+B.shape.locIJ(imin,j);
                c   =Cdata+C.shape.locIJ(imin,j);

                imax=j-1;
                Aend=Adata+  shape.locIJ(imax,j)+1;
                Bend=Bdata+B.shape.locIJ(imax,j)+1;
                Cend=Cdata+C.shape.locIJ(imax,j)+1;

                rangeOp(Op.kind,Op.term,Op.inPlace,C.conjugated!=conjugated,C.conjugated!=B.conjugated,
                        a,Aend,AincT,b,Bend,BincT,c,Cend,CincT,as,bs,cs);

            }
        }
    }
}

UseMatrix UseMatrix::reband(unsigned int SubDiagonal, unsigned int SuperDiagonal, bool Crop) const{
    if(subD()==SubDiagonal and superD()==SuperDiagonal)return *this;
    if(not Crop and (subD()>SubDiagonal or superD()>SuperDiagonal)){
        cout<<"subD,    SubDiagonal "<<subD()<<", "<<SubDiagonal<<endl;
        cout<<"superD,SuperDiagonal "<<superD()<<", "<<SuperDiagonal<<endl;
        cout<<"shape original: "<<strShape()<<endl;
        ABORT("original exceeds new - specify Crop=true for cropping");
    }
    // create new matrix with band width
    symmetry newSym=symm->symm;
    if(SubDiagonal!=SuperDiagonal)newSym=unsymm;
    UseMatrix c(rows(),cols(),SubDiagonal,SuperDiagonal,newSym);
    if(rows()==0 or cols()==0)return c; // empty matrix, no need to reband

    // loop through columns
    complex<double>*ca,*b;
    double *ra;
    for (b=c.cdata;b<c.cdata+c.shape.total;b++)*b=0;

    unsigned int inca=  storeIncI();
    unsigned int incb=c.storeIncI();

    for (unsigned int j=0;j<cols();j++){
        int imin=min(int(rows())-1,max(0,int(j)-int(min(superD(),SuperDiagonal))));
        int imax=min(int(rows())-1,int(j)+int(min(subD(),SubDiagonal)));
        if(cdata!=0){
            for(ca=cdata+shape.locIJ(imin,j),b=c.cdata+c.shape.locIJ(imin,j);ca<cdata+shape.locIJ(imax,j)+1;ca+=inca,b+=incb)*b=*ca;
        } else if(rdata!=0) {
            for(ra=rdata+shape.locIJ(imin,j),b=c.cdata+c.shape.locIJ(imin,j);ra<rdata+shape.locIJ(imax,j)+1;ra++,b++)*b=*ra;
        } else ABORT("cannot reband empty matrix");
    }
    return c;
}

// apply operation to pointer ranges
template<class Ad,class Bd,class Cd,class As,class Bs,class Cs>
void UseMatrix::rangeOp(operation Kind, terms Term, bool InPlace, bool Aconjg, bool Bconjg,
                        Ad a,Ad Aend,unsigned int Ainc,
                        Bd b,Bd Bend,unsigned int Binc,
                        Cd c,Cd Cend,unsigned int Cinc,
                        As as, Bs bs, Cs cs) const
{
    Ad atemp=0;
    Bd btemp=0;

    int ainc,binc;

    if(InPlace){
        // in-place operations
        switch(Term){
        case scalar:
            // scalar operations
            switch(Kind){
            case      add:for(;c<Cend;c+=Cinc)*c+=cs;break;
            case subtract:for(;c<Cend;c+=Cinc)*c-=cs;break;
            case  product:for(;c<Cend;c+=Cinc)*c*=cs;break;
            case quotient:for(;c<Cend;c+=Cinc)*c/=cs;break;
            case   assign:for(;c<Cend;c+=Cinc)*c =cs;break;
            default: ABORT("operation not implemented in-Place scalar");
            }
            break;

        case binary:
            // binary operations


            // get conjugated copy of the data if needed
            binc=Binc;
            if(Bconjg)conjg(b,Bend,binc,btemp);

            switch(Kind){
            case      add:for(;c<Cend and b<Bend;c+=Cinc,b+=binc)*c+=*b;break;
            case subtract:for(;c<Cend and b<Bend;c+=Cinc,b+=binc)*c-=*b;break;
            case  product:for(;c<Cend and b<Bend;c+=Cinc,b+=binc)*c*=(*b);break;
            case quotient:for(;c<Cend and b<Bend;c+=Cinc,b+=binc)*c/=(*b);break;
            case  assign: for(;c<Cend and b<Bend;c+=Cinc,b+=binc)*c=*b;break;
            default: ABORT("operation not implemented inPlace binary");
            }
            // remove storage used for conjugation
            if(Bconjg)delete[] btemp;
            break;
        default:
            ABORT("in-place not implemented");
            cout<<"shut up the compiler about unused variables"<<as<<bs;
        }
    } else {
        switch(Term){
        case ternary:
            // ternary operations

            // get conjugated copy of the data if needed
            ainc=Ainc;
            binc=Binc;
            if(Aconjg)conjg(a,Aend,ainc,atemp);
            if(Bconjg)conjg(b,Bend,binc,btemp);

            switch(Kind){
            case add:      for(;c<Cend;a+=ainc,b+=binc,c+=Cinc)*c=*a+*b;break;
            case subtract: for(;c<Cend;a+=ainc,b+=binc,c+=Cinc)*c=*a-*b;break;
            case product:  for(;c<Cend;a+=ainc,b+=binc,c+=Cinc)*c=*a*(*b);break;
            case quotient: for(;c<Cend;a+=ainc,b+=binc,c+=Cinc)*c=*a/(*b);break;
            default: ABORT("operation not defined for ternary operator");
            }
            // remove storage used for conjugation
            if(Aconjg)delete[] atemp;
            if(Bconjg)delete[] btemp;
            break;

        case scalar:
            // scalar operations

            // get conjugated copy of the data if needed
            binc=Binc;
            if(Bconjg)conjg(b,Bend,binc,btemp);

            switch(Kind){
            case      add:for(;c<Cend;b+=binc,c+=Cinc)*c=*b+cs;break;
            case subtract:for(;c<Cend;b+=binc,c+=Cinc)*c=*b-cs;break;
            case  product:for(;c<Cend;b+=binc,c+=Cinc)*c=*b*cs;break;
            case quotient:for(;c<Cend;b+=binc,c+=Cinc)*c=*b*cs;break;
            default: ABORT("operation not implemented out-of-place scalar");
            }
            if(Bconjg)delete[] btemp;
            break;
        default: ABORT("out-of-place not implemented");
        }
    }
}

// public interface of the operators
UseMatrix   UseMatrix::operator* (const UseMatrix & B) const {UseMatrix C(rows(),B.cols());multiply(B,C,1.,0.);return C;}

// for transition from Eigen usage
Eigen::VectorXcd UseMatrix::operator*(const Eigen::VectorXcd & V){
    if (shape.ncols!=V.rows()) ABORT("cannot multiply, dimensions do not match");
    Eigen::VectorXcd mV;
    mV=Eigen::Map<Eigen::MatrixXcd>(cdata,shape.nrows,shape.ncols)*V;
    return mV;
}

UseMatrix UseMatrix::Constant(unsigned int Nrows, unsigned int Ncols, const std::complex<double> C){
    UseMatrix I(Nrows,Ncols);
    for(unsigned int n=0;n<Nrows*Ncols;n++)*(I.cdata+n)=C;
    return I;
}
UseMatrix UseMatrix::Random(unsigned int Nrows, unsigned int Ncols){
    UseMatrix I(Nrows,Ncols);
    for(unsigned int n=0;n<Nrows*Ncols;n++)*(I.cdata+n)=complex<double>((rand() % 10000)*1.e-4,(rand() % 10000)*1.e-4);
    return I;
}
UseMatrix UseMatrix::RandomReal(unsigned int Nrows, unsigned int Ncols){
    UseMatrix I(Nrows,Ncols);
    for(unsigned int n=0;n<Nrows*Ncols;n++)*(I.cdata+n)=complex<double>((rand() % 10000)*1.e-4,0.);
    return I;
}


//template<class T>
//UseMatrix UseMatrix::FromVector(const T &vector)
//{
//    UseMatrix I(vector.size(), 1);
//    for (unsigned int i=0; i!=vector.size(); ++i) { *(I.cdata+i)=vector.at(i); }
//    return I;
//}

UseMatrix UseMatrix::FromVector(const std::vector<complex<double> > &vector)
{
    UseMatrix I(vector.size(), 1);
    for (unsigned int i=0; i!=vector.size(); ++i) { *(I.cdata+i)=vector.at(i); }
    return I;
}

UseMatrix UseMatrix::FromVector(const std::vector<double> &vector)
{
    UseMatrix I(vector.size(), 1);
    for (unsigned int i=0; i!=vector.size(); ++i) { *(I.cdata+i)=vector.at(i); }
    return I;
}

UseMatrix UseMatrix::Zero(unsigned int Nrows, unsigned int Ncols){
    UseMatrix Z(Nrows,Ncols);
    Z=Constant(Nrows,Ncols,0.);
    return Z;
}
UseMatrix UseMatrix::Identity(unsigned int nrows, unsigned int ncols){
    UseMatrix I(nrows,ncols);
    I=Zero(nrows,ncols);
    for (unsigned int n=0;n<std::min(nrows,ncols)*ncols;n+=ncols+1) *(I.cdata+n)=1.;
    return I;
}

UseMatrix UseMatrix::col(unsigned int I) const {return block(0,I,rows(),1);}
UseMatrix UseMatrix::row(unsigned int I) const {return block(I,0,1,cols());}
UseMatrix UseMatrix::leftCols( unsigned int N){return block(0,       0,rows(),N);}
UseMatrix UseMatrix::rightCols(unsigned int N){return block(0,cols()-N,rows(),N);}
UseMatrix UseMatrix::block(unsigned int I, unsigned int J, unsigned int M, unsigned int N) const {
    if(I+M>rows() or J+N>cols())ABORT("sub-block exceeds matrix, I,J,M,N, dimensions: "
                                      +tools::str(I)+","+tools::str(J)+","+tools::str(M)+","
                                      +tools::str(N)+","+tools::str(rows())+","+tools::str(cols()));
    UseMatrix b;
    b.shape=Shape(M,N,shape.locationIJ,shape.leadDim,shape.total,M,N,I,J,&shape);
    b.shape.locationIJ=any_block;
    b.conjugated=conjugated;
    b.symm->symm=unknown;
    b.rdata=rdata;
    b.cdata=cdata;
    return b;
}

UseMatrix UseMatrix::asDiagonal() const {
    ABORT("this does not work yet");
    if(rows()>1 and cols()>1)ABORT("only single row or single column can be interpreted as diagonal matrix");
    UseMatrix b;
    unsigned int d=max(rows(),cols());
    //    b.shape=Shape(d,d,0,0);
    //    b.shape.locationIJ=asDiagonal_map;
    //    b.shape=Shape(M,N,shape.locationIJ,shape.leadDim,shape.total,M,N,I,J,&shape);
    b.shape=Shape(d,d,asDiagonal_map,1,shape.total,d,d,0,0,0);
    b.conjugated=conjugated;
    b.symm->symm=symmetric;
    b.rdata=rdata;
    b.cdata=cdata;
    b.shape.print("diagonal");
    return b;
}


// view on the band of a matrix
UseMatrix UseMatrix::band(unsigned int SubD, unsigned int SuperD) const{
    UseMatrix b;
    b.shape=Shape(rows(),cols(),any_block,shape.leadDim,shape.total,SubD+1,SuperD+1,0,0,&shape);
    b.conjugated=conjugated;
    b.rdata=rdata;
    b.cdata=cdata;
    return b;
}

void UseMatrix::triFactor() const {
    if(trf==0)const_cast<UseMatrix*>(this)->trf=new TriFactorLU();
    trf->reFactor(*this);
}

UseMatrix UseMatrix::inverse() const {
    UseMatrix D;
    if(rdata==0 and cdata==0)return D; // nothing to be done for emtpy matrix
    if(rdata!=0)ABORT("inverse at present only for complex matrices");
    if(!isShape(full_normal))ABORT("inverse at present only for full_normal matrices");

    //special case dim=1 (Lapack does not seem to handle this right)
    if(rows()==1){
        D=*this;
        D(0,0)=1./D(0,0).complex();
        return D;
    }
    triFactor(); // update the factorization, if needed
    D=trf->inverse();
    return D;
}
UseMatrix & UseMatrix::solve(UseMatrix &Rhs) const {
    if(cols()!=Rhs.rows())ABORT("cannot solve - dimensions do not match: "+tools::str(cols())+"!="+tools::str(Rhs.rows()));
    if(rdata==0 and cdata==0)return Rhs; // nothing to be done for emtpy matrix
    if(rdata!=0)ABORT("inverse at present only for complex matrices");
    if(!isShape(full_normal))ABORT("system solving at present only for full_normal matrices");

    triFactor(); // update the factorization, if needed
    return trf->solve(lapackTrans(),Rhs);
}

char UseMatrix::lapackTrans() const {
    if(shape.storeTranspose()){
        if(conjugated)return 'h';
        else return 't';
    }
    if(conjugated)ABORT("date conjugated, but matrix not transposed, cannot use in lapack");
    return 'n';
}

UseMatrix UseMatrix::transpose() const {
    UseMatrix b;
    b.shape=Shape(cols(),rows(),any_transpose,leadDim(),shape.total,superD()+1,subD()+1,0,0,&shape);
    b.conjugated=conjugated;
    b.rdata=rdata;
    b.cdata=cdata;
    return b;
}
UseMatrix UseMatrix::adjoint() const {
    UseMatrix b;
    b.shape=Shape(cols(),rows(),any_transpose,leadDim(),shape.total,superD()+1,subD()+1,0,0,&shape);
    b.conjugated=not conjugated;
    b.rdata=rdata;
    b.cdata=cdata;
    return b;
}
complex<double> UseMatrix::dot(const UseMatrix &B) const{
    if(shape.nrows!=B.shape.nrows)ABORT("cannot form dot product: dimensions do not match");
    if(shape.ncols!=1 or B.shape.ncols!=1)ABORT("dot product is between column vectors only");
    complex<double> c=0;
    complex<double>*a=cdata,*b=B.cdata;
    a+=  shape.locIJ(0,0);
    b+=B.shape.locIJ(0,0);
    complex<double>*aEnd=a+rows();
    for(;a<aEnd;a++,b++)c+=*a*(*b);
    return c;
}

bool UseMatrix::isIdentity(double eps) const {
    if(shape.nrows!=shape.ncols)return false;
    for(unsigned int j=0;j<cols();j++){
        for(unsigned int i=0;i<j;i++)       if(abs(*(cdata+shape.locIJ(i,j)))>eps)return false;
        for(unsigned int i=j+1;i<rows();i++)if(abs(*(cdata+shape.locIJ(i,j)))>eps)return false;
        if(abs(*(cdata+shape.locIJ(j,j))-1.)>eps)return false;
    }
    return true;
}
bool UseMatrix::isDiagonal(double eps) const {
    if(shape.nrows!=shape.ncols)return false;
    for(unsigned int j=0;j<cols();j++){
        for(unsigned int i=0;i<j;i++)       if(abs(*(cdata+shape.locIJ(i,j)))>eps)return false;
        for(unsigned int i=j+1;i<rows();i++)if(abs(*(cdata+shape.locIJ(i,j)))>eps)return false;
    }
    return true;
}

bool UseMatrix::isZero(double eps) const {
    for(unsigned int j=0;j<cols();j++){
        if       (cdata!=0){
            for(complex<double>* a=begCol(cdata,j);a<endCol(cdata,j);a++)if(abs(a->real())>eps or abs(a->imag())>eps)return false;
        } else if(rdata!=0){
            for(        double*  a=begCol(rdata,j);a<endCol(rdata,j);a++)if(abs(*(a))>eps)return false;
        }
    }
    return true; // includes size zero matrix
}

unsigned int UseMatrix::nonZeros(double Eps) const {
    unsigned int cnt=0;
    for(unsigned int j=0;j<cols();j++){
        if       (cdata!=0){
            for(complex<double>* a=begCol(cdata,j);a<endCol(cdata,j);a++)if(abs(a->real()) or abs(a->imag())>Eps)cnt++;
        } else if(rdata!=0){
            for(        double*  a=begCol(rdata,j);a<endCol(rdata,j);a++)if(abs(*(a))>Eps)cnt++;
        }
    }
    // duplicate for testing
    unsigned int nonZ,trueSub,trueSuper;
    char kind;
    diagnose(Eps,nonZ,trueSub,trueSuper,kind);
    if(cnt!=nonZ){
        print("matrix",2);
        ABORT("diagnose and nonZeros differ");
    }
    return cnt;
}

vector<unsigned int> UseMatrix::locNonZero(double Eps) const{
    vector<unsigned int> loc(2);
    for (loc[1]=0;loc[1]<cols();loc[1]++)
        for (loc[0]=0;loc[0]<rows();loc[0]++)
            if(abs(cdata[shape.locIJ(loc[0],loc[1])])>Eps)return loc;
    return loc;
}


// row-wise and column-wise maximal norm
double UseMatrix::maxNorm(std::vector<double>& RowMax, std::vector<double>& ColMax) const {

    double maxNrm=0.;
    ColMax.resize(cols());
    RowMax.assign(rows(),0.);
    for(unsigned int j=0;j<cols();j++){
        double colMax=0.;
        unsigned int aBeg=shape.begCol(j),aEnd=shape.endCol(j);
        unsigned int i=j+aBeg-shape.locIJ(j,j);
        if(cdata!=0)
            for(const complex<double>* a=cdata+aBeg;a<cdata+aEnd;a++,i++){
                double nrm=norm(*a);
                if(nrm>colMax   )colMax   =nrm;
                if(nrm>RowMax[i])RowMax[i]=nrm;
            }
        else
            for(const double* a=rdata+aBeg;a<rdata+aEnd;a++,i++){
                double nrm=abs(*a);
                if(nrm>colMax   )colMax   =nrm;
                if(nrm>RowMax[i])RowMax[i]=nrm;
            }

        ColMax[j]=colMax;
        maxNrm=max(colMax,maxNrm);
    }
    return maxNrm;
}

// row-wise and column-wise maximal norm
double UseMatrix::maxRealImag(std::vector<double>& RowMax, std::vector<double>& ColMax) const {

    double maxNrm=0.;
    ColMax.resize(cols());
    RowMax.assign(rows(),0.);
    for(unsigned int j=0;j<cols();j++){
        double colMax=0.;
        unsigned int aBeg=shape.begCol(j),aEnd=shape.endCol(j);
        unsigned int i=j+aBeg-shape.locIJ(j,j);
        if(cdata!=0)
            for(const complex<double>* a=cdata+aBeg;a<cdata+aEnd;a++,i++){
                double nrm=max(abs(a->real()),abs(a->imag()));
                if(nrm>colMax   )colMax   =nrm;
                if(i<0 or i>=RowMax.size())
                    Sstr+"shape"+shape.str()+i+j+rows()+cols()+(aEnd-aBeg)+Sendl;
                if(nrm>RowMax[i])RowMax[i]=nrm;
            }
        else
            for(const double* a=rdata+aBeg;a<rdata+aEnd;a++,i++){
                double nrm=abs(*a);
                if(nrm>colMax   )colMax   =nrm;
                if(nrm>RowMax[i])RowMax[i]=nrm;
            }

        ColMax[j]=colMax;
        maxNrm=max(colMax,maxNrm);
    }
    return maxNrm;
}

/// RealMat storage of same shape as in complex matrix
/// returns false and zero size Phase and RealMat if decomposition not possible
bool UseMatrix::polar(std::vector<std::complex<double> > &Phase, std::vector<double> &RealMat) const{
    Phase.assign(rows(),1.e-12);
    RealMat.resize(shape.total);
    for(unsigned int j=0;j<cols();j++){
        // get largest element in column
        unsigned int i=j+shape.begCol(j)-shape.locIJ(j,j);
        for(complex<double> *a=begCol(cdata,j);a<endCol(cdata,j);a++,i++)
            if(abs(Phase[i].real())<abs(a->real()) and abs(Phase[i].imag())<abs(a->imag()))Phase[i]=*a;
    }
    vector<double>maxRow;
    for(unsigned int i=0;i<Phase.size();i++){
        maxRow.push_back(abs(Phase[i]));
        Phase[i]/=max(1.e-12,maxRow.back());
    }

    for(unsigned int j=0;j<cols();j++){
        // check phases in column
        unsigned int i=j+shape.begCol(j)-shape.locIJ(j,j);
        for(complex<double> *a=begCol(cdata,j);a<endCol(cdata,j);a++,i++){
            complex<double> c=*a*conj(Phase[i]);
            if(pow(c.imag(),2)<=norm(c)*1e-24)
                RealMat[shape.locIJ(i,j)]=c.real();
            else {
                //                cout<<"failed "<<i<<" "<<c<<endl;
                // does not have polar structure, terminate:
                Phase.clear();
                RealMat.clear();
                return false;
            }
        }
    }
    return true;
}


// fast diagnosis of matrix structure
void UseMatrix::diagnose(double Eps, unsigned int &NonZero, unsigned int &TrueSub, unsigned int &TrueSuper, char &DataType) const {

    if(size()==0){
        NonZero=0;
        TrueSub=0;
        TrueSuper=0;
        DataType='c';
        return;
    }

    vector<double>epsI(rows(),0.),epsJ(cols(),0.);
    if(Eps!=0.){
        maxRealImag(epsI,epsJ);
        for(unsigned int k=0;k<rows();k++)epsI[k]*=Eps;
        for(unsigned int k=0;k<cols();k++)epsJ[k]*=Eps;
    }

    NonZero=0;
    TrueSub=0;
    TrueSuper=0;
    if(cdata!=0){
        if(     cdata->imag()<=max(epsI[0],epsJ[0]))DataType='r';
        else if(cdata->real()<=max(epsI[0],epsJ[0]))DataType='i';
        else DataType='c';

        for(unsigned int j=0;j<cols();j++){
            unsigned int iBeg=shape.begCol(j);
            unsigned int i=j+iBeg-shape.locIJ(j,j);
            const complex<double> *aNonz=0,*aZero=cdata+iBeg-1,*aDiag=cdata+shape.locIJ(j,j);
            for(const complex<double>* a=cdata+iBeg;a<cdata+shape.endCol(j);a++,i++){
                if(std::isnan(a->real()) or std::isnan(a->imag())){
                    print("not a number",2);
                    ABORT("matrix contains non-numbers");
                }
                double epsij=max(epsJ[j],epsI[i]);
                if((abs(a->real())>epsij or abs(a->imag())>epsij))
                {


                    NonZero++;
                    aNonz=a;
                    // update data type (unless already complex)
                    if(DataType!='c'){
                        if(DataType=='r'){
                            if(a->imag()>epsij)DataType='c';
                        } else
                            if(a->real()>epsij)DataType='c';
                    }
                }
                else
                    if(aNonz==0)aZero=a;
            }

            // update lower and upper band width from aZero and aNonz
            if(aNonz==0)aNonz=aDiag;
            TrueSuper=max(int(TrueSuper),int((aDiag-aZero)-1));
            TrueSub  =max(int(TrueSub),  int((aNonz-aDiag))  );
        }
    }

    else {
        DataType='r';
        for(unsigned int j=0;j<cols();j++){
            unsigned int iBeg=shape.begCol(j);
            unsigned int i=j+iBeg-shape.locIJ(j,j);
            const double *aNonz=0,*aZero=rdata+iBeg-1,*aDiag=rdata+shape.locIJ(j,j);
            for(const double* a=rdata+iBeg;a<rdata+shape.endCol(j);a++,i++){
                double epsij=max(epsJ[j],epsI[i]);
                if(abs(*a)>epsij){
                    NonZero++;
                    aNonz=a;
                }
                else
                    if(aNonz==0)aZero=a;
            }

            // update lower and upper band width from aZero and aNonz
            if(aNonz==0)aNonz=aDiag;
            TrueSuper=max(int(TrueSuper),int((aDiag-aZero)-1));
            TrueSub  =max(int(TrueSub),  int((aNonz-aDiag))  );
        }
    }
}

double UseMatrix::maxAbsVal() const {
    double maxAbs=0.;
    if       (cdata!=0){
        for(unsigned int j=0;j<cols();j++)
            for(complex<double>* a=begCol(cdata,j);a<endCol(cdata,j);a++)
                maxAbs=max(maxAbs,norm(*a));
    } else if(rdata!=0){
        for(unsigned int j=0;j<cols();j++)
            for(       double*  a=begCol(rdata,j);a<endCol(rdata,j);a++)
                maxAbs=max(maxAbs,abs(*a));
    }
    if(cdata!=0)return sqrt(maxAbs);
    else        return maxAbs;
}

vector<complex<double> > UseMatrix::extractSubmatrix(std::vector<unsigned int> Rows, std::vector<unsigned int> Cols) const{
    vector<complex<double> > dat(Rows.size()*Cols.size());
    for(unsigned int j=0,k=0;j<Cols.size();j++)
        for(unsigned int i=0;i<Rows.size();i++,k++){
            dat[k]=*(cdata+shape.locIJ(Rows[i],Cols[j]));
        }
    return dat;
}

bool UseMatrix::operator==(const UseMatrix &rhs) const {return  shape==rhs.shape and compare(rhs,0.);}

bool UseMatrix::compare(const UseMatrix &rhs,double Eps) const {
    if(this->cols() != rhs.cols() or this->rows() != rhs.rows()) return false;
    for(unsigned int j=0;j<this->cols();j++){
        double epsAbs=Eps;
        if(Eps!=0)epsAbs=Eps*std::max(col(j).maxAbsVal(),rhs.col(j).maxAbsVal());
        if(this->cdata!=0){
            if(rhs.cdata==0) return false;
            complex<double>* b=rhs.begCol(rhs.cdata,j);
            for(complex<double>* a=this->begCol(this->cdata,j);a<this->endCol(this->cdata,j);a++,b++){
                if(b>=rhs.endCol(rhs.cdata,j)) return false;
                if(abs(a->real()-b->real())+abs(a->imag()-b->imag())>epsAbs)  return false;
            }
        } else {
            ABORT("UseMatrix::operator== real not implemented");
        }
    }
    return true;
}

void UseMatrix::show(string text) const {print(text,0);}

void UseMatrix::print(string text, unsigned int Digits) const {
    cout<<str(text,Digits)<<endl;
}

string UseMatrix::str(string text, unsigned int Digits) const {
    ostringstream oss;
    // column numbering
    if(text==""); // do not print header
    else if(text.length()>9)oss<<text<<endl<<"row\\col  ";
    else oss<<"\n"<<setw(text.length())<<text<<setw(9-text.length())<<"";
    if(Digits>0)
        for (unsigned int n=0;n<cols();n++){
            if(n==0)oss<<setw(Digits+1)<<0;
            else oss<<setw(Digits+7)<<n;
        }
    oss<<endl;

    for (unsigned int m=0;m<rows();m++){
        // row numbering
        oss<<setw(4)<<m;
        // real parts
        for (unsigned int n=0;n<cols();n++){
            if(Digits>0)oss<<setprecision(Digits)<<setw(Digits+7)<<(*this)(m,n).real();
            else        oss<<tools::zero((*this)(m,n));
        }
        oss<<endl;

        // imaginary parts
        if(Digits>0 and cdata!=0){
            oss<<setw(4)<<" ";
            if(Digits>0)
                for (unsigned int n=0;n<cols();n++)
                    oss<<setprecision(Digits)<<setw(Digits+7)<<(*this)(m,n).imag();
            oss<<endl;
        }
    }
    oss<<shape.str();
    return oss.str();
}

bool UseMatrix::zeroData(const unsigned int Loc, const double Eps) const {
    if     (cdata!=0)return abs((cdata+Loc)->real())<=Eps and abs((cdata+Loc)->imag())<=Eps ;
    else if(rdata!=0)return abs(*(rdata+Loc))<=Eps;
    else return true;
}
bool UseMatrix::realData(const unsigned int Loc, const double Eps) const {
    if(cdata!=0) return abs(imag(*(cdata+Loc)))<=Eps;
    else return true;

}
bool UseMatrix::equalData(const unsigned int LocA, const unsigned int LocB, const double Eps) const {
    if     (cdata!=0)return abs((cdata+LocA)->real()-(cdata+LocB)->real())<=Eps and abs((cdata+LocA)->imag()-(cdata+LocB)->imag())<=Eps;
    else if(rdata!=0)return abs(*(rdata+LocA)-*(rdata+LocB))<=Eps;
    else return true;
}
bool UseMatrix::conjgData(const unsigned int LocA, const unsigned int LocB, const double Eps) const {
    if     (cdata!=0)return abs((cdata+LocA)->real()-(cdata+LocB)->real())<=Eps and abs((cdata+LocA)->imag()+(cdata+LocB)->imag())<=Eps;
    else if(rdata!=0)return abs(*(rdata+LocA)-     *(rdata+LocB) )<=Eps;
    else return true;
}

unsigned int UseMatrix::trueSubD(double Eps) const {
    unsigned int nonZ,trueSub,trueSuper;
    char kind;
    diagnose(Eps,nonZ,trueSub,trueSuper,kind);
    return trueSub;
}

unsigned int UseMatrix::trueSuperD(double Eps) const {
    // duplicate for testing
    unsigned int nonZ,trueSub,trueSuper;
    char kind;
    diagnose(Eps,nonZ,trueSub,trueSuper,kind);
    return trueSuper;
}

bool UseMatrix::isHermitian(const double Eps, double Threshold) const {
    if(symm->symm==hermitian)return true;
    if(symm->symm==symmetric or symm->symm==unsymm)return false;
    if(rows()!=cols())return false;
    for(unsigned int j=0;j<cols();j++){
        double nrm=col(j).maxAbsVal();
        if(nrm<=Threshold and row(j).maxAbsVal()<=Threshold)continue;
        double epsMax=max(Eps*col(j).maxAbsVal(),Threshold);
        for(unsigned int i=j;i<min(rows(),min(j+subD()+1,j+superD()+1));i++){
            if(not conjgData(shape.locIJ(i,j),shape.locIJ(j,i),epsMax)){
                return false;
            }
        }
    }
    symm->symm=hermitian;
    return true;
}
bool UseMatrix::isSymmetric( const double Eps) const {
    // allowed deviation from symmetry relative to the maximal value in the column
    if(symm->symm==symmetric)return true;
    if(symm->symm==hermitian or symm->symm==unsymm)return false;
    if(rows()!=cols())return false;
    // check
    for(unsigned int j=0;j<cols();j++){
        double epsMax=Eps*(this->col(j).maxAbsVal());
        for(unsigned int i=j;i<min(rows(),min(j+subD()+1,j+superD()+1));i++){
            if(not equalData(shape.locIJ(i,j),shape.locIJ(j,i),epsMax)) {
                return false;
            }
        }
    }
    symm->symm=symmetric; // remember for later
    return true;
}
bool UseMatrix::isReal( const double Eps) const {
    for(unsigned int j=0;j<cols();j++) {
        double epsMax=Eps*(this->col(j).maxAbsVal());
        if(j<rows())epsMax=max(epsMax,Eps*this->row(j).maxAbsVal());
        for(unsigned int i=0;i<rows();i++) {
            if(not realData(shape.locIJ(i,j),epsMax)) return false;
        }
    }
    return true;
}

string UseMatrix::type(double BandRatio) const {
    string t;

    unsigned int nonZ,subd,superd;
    char data;
    diagnose(1.e-14,nonZ,subd,superd,data);

    switch (data){
    case 'r': t+="r";break;
    case 'i': t+="i";break;
    case 'c': t+="c";break;
    }

    // hermiticity/symmetry
    if(     isSymmetric())t+="s";
    else if(isHermitian())t+="h";
    else                  t+="g";

    if(nonZ<min(rows(),cols()) and isDiagonal(0.))
        t+="d";
    else if(double(nonZ)<=max(BandRatio*double(rows()*cols()),double(min(rows(),cols())))){
        // diagonal or banded (defaults to full)
        if (rows()==cols()){
            if(subd==0 and superd==0)
                t+="d";
            else if(subd+1<rows() and superd+1<cols()){
                if(subd==0 and superd==0)//debug: keep for now, although superseeded by above
                    t+="d";
                else if(double(subd+superd+1)*min(rows(),cols())<=BandRatio*double(rows()*cols()))
                    t+="b";
            }
        }
        // general sparsity
        else
            t+="s";
    }

    return t;
}

vector<complex<double> > UseMatrix::extractDiagonal(){
    vector<complex<double> > res(min(rows(),cols()));
    for(unsigned int k=0;k<res.size();k++)res[k]=*(cdata+shape.locIJ(k,k));
    return res;
}
vector<complex<double> > UseMatrix::extractSubmatrix(const std::vector<unsigned int> &Rows, const std::vector<unsigned int> &Cols){
    vector<complex<double> > res(Rows.size()*Cols.size());
    for(unsigned int j=0,ij=0;j<Cols.size();j++)
        for(unsigned int i=0;i<Rows.size();i++,ij++)
            res[ij]=*(cdata+shape.locIJ(i,j));
    return res;
}

UseMatrix &UseMatrix::purge(double Eps,double EpsAbs){

    if(Eps<=0. and EpsAbs<=0.)return *this; // zero threshold, no purge

    if(not isBasic() and not shape.storeFull())ABORT("cannot .purge(...) non-basic matrices, is: "+strShape());

    vector<double> epsI,epsJ;
    if(cdata!=0){
        maxRealImag(epsI,epsJ);
        for(unsigned int k=0;k<epsI.size();k++)epsI[k]=max(epsI[k]*Eps,EpsAbs);
        for(unsigned int k=0;k<epsJ.size();k++)epsJ[k]=max(epsJ[k]*Eps,EpsAbs);
        for(unsigned int j=0;j<cols();j++) {
            unsigned int i=j+shape.begCol(j)-shape.locIJ(j,j);
            for(complex<double>*a=cdata+shape.begCol(j);a<cdata+shape.endCol(j);a++,i++) {
                double epsij=max(epsI[i],epsJ[j]);
                if(abs(a->imag())<epsij)*a=complex<double>(a->real(),0.);
                if(abs(a->real())<epsij)*a=complex<double>(0.,a->imag());
            }
        }
    }

    else {
        maxNorm(epsI,epsJ);
        for(unsigned int k=0;k<epsI.size();k++)epsI[k]=max(sqrt(epsI[k])*Eps,EpsAbs);
        for(unsigned int k=0;k<epsJ.size();k++)epsJ[k]=max(sqrt(epsJ[k])*Eps,EpsAbs);
        for(unsigned int j=0;j<cols();j++) {
            unsigned int i=j+shape.begCol(j)-shape.locIJ(j,j);
            for(double*a=rdata+shape.begCol(j);a<rdata+shape.endCol(j);a++,i++)
                if(abs(*a)<max(epsI[i],epsJ[j]))*a=0;
        }
    }
    return *this;
}

UseMatrix UseMatrix::expandConst() const{return reband(rows()-1,cols()-1);}

UseMatrix & UseMatrix::expand(){
    UseMatrix c(*this);
    c=reband(rows()-1,cols()-1);
    swap(c);
    return *this;
}
UseMatrix & UseMatrix::compress(double Eps, bool Symm){
    STARTDEBUG(compress);
    if(Symm)storeSymmetricBand(Eps);
    else storeBand(Eps);
    STOPDEBUG(compress);
    return *this;
}

unsigned int UseMatrix::storeIncI() const{
    unsigned int i=1;
    if(shape.storeTranspose())i=shape.leadDim;
    if(i>1 and isBand())i--;
    return i;
}

void UseMatrix::storeSymmetricBand(double Eps){
    if(not isSymmetric(Eps))ABORT("not symmetric - cannot store lower band, tolerance = "+tools::str(Eps));
    int subd=trueSubD(Eps);
    // create new
    UseMatrix c(rows(),cols(),subd,subd,symmetric);
    c=this->band(subd,subd);
    swap(c);
}
void UseMatrix::storeBand(double Eps){
    int subd=trueSubD(Eps),superd=trueSuperD(Eps);
    // create new
    UseMatrix c(rows(),cols(),subd,superd);
    c=this->band(subd,superd);
    swap(c);
}


//void UseMatrix::Usage(){
//    using namespace Eigen;

//    UseMatrix m(3,2);
//    for (unsigned int i=0;i<m.rows();i++){
//        for (unsigned int j=0;j<m.cols();j++){
//            if(i<j)m(i,j)=complex<double>(1+i+0.1*(1+j),0.); // this is OK, casting works
//            else m(i,j)=1+i+0.01*(1+j);
//        }
//    }
//    m.print("real");


//    UseMatrix c(3,2);
//    for (unsigned int i=0;i<c.rows();i++){
//        for (unsigned int j=0;j<c.cols();j++){
//            if(i<j)c(i,j)=complex<double>(1+i+0.01*(1+j),-1.*((double) (j+ 0.1*i)));
//            else c(i,j)=1+i+0.01*(1+j);
//        }
//    }
//    c.print("complex",6);

//    c*=2.;
//    c.print("c*2");

//    c*=complex<double>(0.,1.);
//    c.print("i*c");

//    UseMatrix C(c);
//    C.print("copy");

//    cout<<"\nthis works:"<<endl;
//    double d=c(1,0).real();
//    cout<<"value of c(1,0) = "<<d<<endl;

//    UseMatrix b;
//    b=m.block(1,0,2,2);
//    b.print("block");

//    m.block(1,1,1,2).print("also direct printing of derived matrix works");

//    m.transpose().print("...and direct printing of transpose");

//    // check transposition
//    UseMatrix t;
//    UseMatrix tt;
//    m=UseMatrix::Random(17,32);
//    t=m.transpose();
//    tt=t.transpose();
//    if(not (m-tt).isZero())ABORT("FAIL double transpose");
//    cout<<"OK double transpose "<<m.rows()<<" x "<<m.cols()<<endl;

//    // check normal matrix multiply
//    UseMatrix R=UseMatrix::Random(30,44),S=UseMatrix::Random(44,22);
//    UseMatrix RS=R*S;
//    MatrixXcd eRS=Eigen::Map<MatrixXcd>(R.data(),R.rows(),R.cols())*Eigen::Map<MatrixXcd>(S.data(),S.rows(),S.cols());
//    if(not (eRS-Eigen::Map<MatrixXcd>(RS.data(),RS.rows(),RS.cols())).isZero())ABORT("FAIL matrix multiply ");
//    cout<<"OK matrix multiply "<<R.rows()<<" x "<<R.cols()<<" x "<<S.cols()<<endl;

//    // test inverse
//    UseMatrix A=UseMatrix::Random(17,17);
//    UseMatrix invA=A.inverse();
//    if(not (A*invA).isIdentity(1.e-14))ABORT("FAIL inverse");
//    cout<<"OK inverse "<<A.rows()<<" x "<<A.cols()<<endl;

//}
//void UseMatrix::Test(){
//    using namespace Eigen;

//    int const d0=5,d1=17,d2=32,d3=44;
//    UseMatrix m;

//    // check transposition
//    UseMatrix t;
//    UseMatrix tt;
//    m=UseMatrix::Random(d1,d2);
//    t=m.transpose();
//    tt=t.transpose();
//    if(not (m-tt).isZero())ABORT("FAIL double transpose");
//    cout<<"OK double transpose "<<m.rows()<<" x "<<m.cols()<<endl;

//    // check normal matrix multiply
//    UseMatrix R=UseMatrix::Random(d2,d3),S=UseMatrix::Random(d3,d1);
//    UseMatrix RS=R*S;
//    MatrixXcd eRS=Eigen::Map<MatrixXcd>(R.data(),R.rows(),R.cols())*Eigen::Map<MatrixXcd>(S.data(),S.rows(),S.cols());
//    if(not (eRS-Eigen::Map<MatrixXcd>(RS.data(),RS.rows(),RS.cols())).isZero())ABORT("FAIL matrix multiply ");
//    cout<<"OK matrix multiply "<<R.rows()<<" x "<<R.cols()<<" x "<<S.cols()<<endl;

//    // test inverse
//    UseMatrix A=UseMatrix::Random(17,17);
//    UseMatrix invA;
//    invA=A.inverse();
//    if(not (A*invA).isIdentity(1.e-14))ABORT("FAIL inverse");
//    cout<<"OK inverse "<<A.rows()<<" x "<<A.cols()<<endl;

//    // test banded
//    m=UseMatrix::Random(d2,d3);
//    UseMatrix b;
//    b=m.band(d0,d0);
//    b.compress();
//    m=b;
//    m.expand();
//    if(not (b-m).isZero())(b-m).show("error in difference band-full");
//    if(not (m-b).isZero())(m-b).show("error in differenc full-band");
//    cout<<"OK compress/expand band matrix"<<endl;

//    // test band on full matrix multiply
//    UseMatrix bf,mf,f=UseMatrix::Random(d3,d0);
//    mf=m*f;
//    bf=b*f;
//    if(not (mf-bf).isZero())(m-b).show("error in difference full*full-band*full");
//    cout<<"OK band*full matrix"<<endl;

//    // generate a complex symmetric matrix
//    UseMatrix zs=UseMatrix::Random(d2,d2);
//    UseMatrix inv=zs.inverse();
//    if(not((inv*zs).isIdentity(1.e-12)))(inv*zs).print("inverse failed");
//    cout<<"OK full general inverse"<<endl;

//    zs=zs.transpose()*zs;
//    if(not zs.isSymmetric(1.e-12))zs.show("not complex symmetric");
//    cout<<"OK complex symmetric"<<endl;

//}

UseMatrix::UseMatrix(const Eigen::MatrixXcd &mat)
    : rdata(0),cdata(0),conjugated(false),symm(new Symmetry(unknown))
{
    trf=new TriFactorLU();
    shape=Shape(mat.rows(),mat.cols());
    allocateData(true);
    for(size_t k=0;k<size();k++)data()[k] = mat.data()[k];
}

