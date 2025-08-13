// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorTensor.h"
#include "tools.h"
#include <typeinfo>

//resolve forward declarations
#include "qtEigenDense.h"
#include <Core>
#ifdef _HAS_LAPACKE_
#include "lapacke.h"
#endif
#include "index.h"
#include "coefficients.h"
#include "coefficientsFloor.h"
#include "parameters.h"
#include "operatorData.h"
#include "printOutput.h"
#ifdef _CBLAS_
#include "cblas.h"
#endif

using namespace std;
using namespace tools;
#include "eigenNames.h"


TIMER(setupTensor1,)
TIMER(setupTensor2,)
OperatorTensor::OperatorTensor(const string & Definition, const Discretization *IDisc, const Discretization *JDisc,
                               const Index *  iFloor, const Index * jFloor, complex<double> Multiplier):
    OperatorSingle("no name",Definition,iFloor,jFloor),mStruc(Fu){
    STARTDEBUG(setupTensor1);

    string hash=Definition;
    hash+=iFloor->hash()+jFloor->hash()+tools::str(Multiplier);
    vector<UseMatrix> matrices;
    OperatorData::operatorData(Definition,iFloor,jFloor,matrices);

    for (unsigned int k=0;k<matrices.size();k++){
        if(k==0)matrices[0]*=Multiplier;
        if(matrices[k].isZero(1e-12)){
            mats.resize(0);
            break;
        }

        // compress to banded storage
        if(matrices.size()==1
                and matrices[k].rows()+matrices[k].cols()>5){
            matrices[k].compress(1.e-12,false);
            if(not matrices[k].isBandNormal())matrices[k].expand();
        }
        const UseMatrix* temp = OperatorSingle::matsAdd(matrices[k],hash+tools::str(k));
        mats.push_back(temp);
    }
    setStructureFlags(Definition);

    STOPDEBUG(setupTensor1);
}

// construct from vector of tensor factors const
OperatorTensor::OperatorTensor(const string & Definition, const Index *IIndex, const Index * JIndex, vector<UseMatrix> & matrices) :
    OperatorSingle("no name",Definition,IIndex,JIndex){

    STARTDEBUG(setupTensor2);
    // check for zero matrices, if non-zero add to the matrix-repository
    string hash=Definition;
    hash+=IIndex->hash()+JIndex->hash();
    for (unsigned int k=0;k<matrices.size();k++){
        if(matrices[k].isZero(1e-14)){
            mats.resize(0);
            break;
        }
        const UseMatrix* temp = OperatorSingle::matsAdd(matrices[k],hash+tools::str(k));
        mats.push_back(temp);
    }

    setStructureFlags(Definition);
    setNorm();
    STOPDEBUG(setupTensor2);
}

void OperatorTensor::setStructureFlags(const string &Definition)
{
    // set structure flags
    if(mats.size()==1){
        mStruc=Fu;
        if(     mats[0]->isIdentity())mStruc=Id;
        else if(mats[0]->isBandNormal())mStruc=Bd;
        if(mats.size()>1 and mats[0]->isBand()!=(mStruc==Bd))ABORT("really bad: "+Definition);
    }
    else if(mats.size()==2){
        mStruc=FuFu;
        if(mats[0]->isIdentity()){
            if(mats[1]->isIdentity())mStruc=Id;
            else mStruc=IdFu;
        } else if(mats[1]->isIdentity())
            mStruc=FuId;
    }
}

void OperatorTensor::matrix(UseMatrix &Mat) const{
    if(mats.size()==1){
        Mat=*mats[0];
    } else
        OperatorSingle::matrix(Mat);
}

// replace all tensor factors by their inverses
void OperatorTensor::inverse(){

    if (definition == "highEVket" or definition == "highEVbra")
        ABORT("OperatorTensor::inverse(): never call for high ev projectors");

    definition += "^-1";
    string hash=definition;
    hash+=iIndex->hash()+jIndex->hash();
    for (unsigned int k = 0; k < mats.size(); k++){
        if (mats[k]->isZero(1.e-14)){
            cout << definition << " is zero ... " << k << endl;
            ABORT("OperatorTensor::inverse");
        }

        //NOTE: eventually, UseMatrix should be tought to handle this without expand
        UseMatrix temp1(mats[k]->expandConst());
        temp1=temp1.inverse();
        mats[k] = OperatorSingle::matsAdd(temp1,hash+tools::str(k));
        if(mats.size()==1)mStruc=Fu;
        if(mats.size()==2)mStruc=FuFu;
    }
}

TIMERSAMPLE(apply,)
void OperatorTensor::apply(std::complex<double> *InOut) const{
    vector<complex<double> >c(InOut,InOut+iIndex->sizeStored());
    apply(c);
    for(unsigned int k=0;k<c.size();k++)InOut[k]=c[k];
}

void OperatorTensor::apply(vector<complex<double> > & InOut) const //Inout = Operator*InOut;
{
    STARTDEBUG(apply);
    switch (mStruc) {
    case Id:
    {
        STOPDEBUG(apply);
        return;
    }
    case Bd:
    {
#ifdef _CBLAS_
        complex<double> zero=0.,one=1.;
        // lapack complex banded by vector
        cout<<"mats[0] "<<mats[0]->strShape()<<" "<<mats[0]->subD()<<" "<<mats[0]->subD()<<" "<<mats[0]->leadDim()<<endl;
        ABORT("case not implemented");
        cblas_zgbmv(CblasColMajor,mats[0]->cblasTrans(),mats[0]->rows(),mats[0]->cols(),mats[0]->subD(),mats[0]->superD(),
                reinterpret_cast<double*>(&one),reinterpret_cast<double*>(mats[0]->data()),mats[0]->leadDim(),
                reinterpret_cast<double*>(InOut.data()),1,reinterpret_cast<double*>(&zero),reinterpret_cast<double*>(InOut.data()),1);
        ABORT("this may not be working");
#else
        DEVABORT("no banded multiply available");
#endif
        break;
    }
    case Fu:
    {
        if(mats.size()==2)ABORT("wrong dimension");
        (Map<VectorXcd>(InOut.data(),InOut.size())) =
                //CHANGE                = (*factor)
                1.
                * Map<MatrixXcd>(mats[0]->data(),mats[0]->rows(),mats[0]->cols())
                * Map<VectorXcd>(InOut.data(),InOut.size());
        break;
    }
    case FuFu:
    {
        Map<MatrixXcd>(InOut.data(),mats[1]->rows(),mats[0]->rows())
                //CHANGE                = (*factor)
                = 1.
                * Map<MatrixXcd>(mats[1]->data(),mats[1]->rows(),mats[1]->cols())
                * Map<const MatrixXcd>(InOut.data(),mats[1]->cols(), mats[0]->cols())
                * (Map<MatrixXcd>(mats[0]->data(),mats[0]->rows(), mats[0]->cols())).transpose();
        break;
    }
    case IdFu:
    {
        Map<MatrixXcd>(InOut.data(),mats[1]->rows(),mats[0]->rows())
                //CHANGE                = (*factor)
                = 1.
                * Map<MatrixXcd>(mats[1]->data(),mats[1]->rows(),mats[1]->cols())
                * Map<const MatrixXcd>(InOut.data(),mats[1]->cols(), mats[0]->cols());
        break;
    }
    case FuId:
    {
        Map<MatrixXcd>(InOut.data(),mats[1]->rows(),mats[0]->rows())
                //CHANGE                = (*factor)
                = 1.
                * Map<const MatrixXcd>(InOut.data(),mats[1]->cols(), mats[0]->cols())
                * (Map<MatrixXcd>(mats[0]->data(),mats[0]->rows(), mats[0]->cols())).transpose();
        break;
    }
    default:
        ABORT("TensorOperator apply not implemented");
    }
    STOPDEBUG(apply)
}

bool OperatorTensor::absorb(OperatorSingle *&Other){
    // only operators with equal indices and equal tensor can be absorbed
    if(iIndex!=Other->iIndex)return false;
    if(jIndex!=Other->jIndex)return false;
    //CHANGE    if(factor!=Other->factor)return false;
    if(typeid(*Other)!=typeid(*this))return false; // can only absorb equal types

    if(mats.size()!=Other->mats.size())return false;
    if(mats.size()>1)return false;

    if(name.find("(fused)")!=string::npos)name+="(fused)";
    definition+=" "+Other->definition;

    UseMatrix* mat=new UseMatrix(mats[0]->rows(),mats[0]->cols());
    *mat=*(mats[0])+*(Other->mats[0]);
    mStruc=Fu;
    mat->expand();
    mat->purge(1.e-12,1.e-12);
    mats[0]=mat;
    return true;
}

string OperatorTensor::str() const {return name+": "+tools::str(mats[0])+" "+tools::str(mats.size());}
string OperatorTensor::strDataStructure() const {
    string s=definition;
    if(mats.size()==0)return s+" (emtpy)";
    for(unsigned int k=0;k<mats.size();k++){
        s+=" "+tools::str(mats[k]->rows())+"x"+tools::str(mats[k]->cols());
    }
    return s;
}

bool OperatorTensor::isZero(double Eps) const {
    /// is zero if any factor is zero or has no factors at all
    for(unsigned int k=0;k<mats.size();k++) {
        if(mats[k]!=0 and mats[k]->isZero(Eps))return true;
    }
    return mats.size()==0;
}

OperatorTensor::~OperatorTensor(){}
