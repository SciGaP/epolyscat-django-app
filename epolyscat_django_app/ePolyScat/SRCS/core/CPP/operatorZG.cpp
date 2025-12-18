// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorZG.h"

#include "qtEigenSparse.h"
//#include <Core>

#include "tools.h"
#include "readInput.h"
#include "useMatrix.h"
#include "index.h"
#include "discretization.h"
#include "operatorSingle.h"
//#include "operator.h"

using namespace std;
#include "eigenNames.h"


OperatorZG::OperatorZG(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf)
    :OperatorFloor("ZG"){
    unpackBasic(Info,Buf);
    dat=addComplex(hashString(_rows,_cols),vector<complex<double> >(Buf.begin(),Buf.begin()+Info[3]));
}

OperatorZG::OperatorZG(const UseMatrix* Dat, string Kind)
    :OperatorFloor(Dat->rows(),Dat->cols(),"ZG")
{

    vector<complex<double> > mDat(Dat->rows()*Dat->cols());

    unsigned int k=0;
    for(unsigned int n=0;n<Dat->cols();n++)
        for(unsigned int m=0;m<Dat->rows();m++,k++)
            mDat[k]=Dat->operator()(m,n).complex();

    string hash=Kind+tools::str(Dat->rows())+tools::str(Dat->cols());
    dat=addComplex(hash,mDat);
    _rows=Dat->rows();
    _cols=Dat->cols();
    oNorm=Dat->maxAbsVal();
}

void OperatorZG::pack(vector<int> & Info,vector<std::complex<double> > &Buf) const {
    Buf.insert(Buf.end(),dat->begin(),dat->end());
    packBasic(Info,Buf);
}

void OperatorZG::axpy(const std::complex<double>&Alfa, const std::complex<double>* X, unsigned int SizX,
                      const std::complex<double>&Beta, std::complex<double>* Y, unsigned int SizY) const{

    // after some consulting of Eigen advice, Y +=Alfa * Matrix * X seems to be the right thing to do
    if(Alfa==0.){
        scale(Beta,Y,SizY);
        return;
    }

    complex<double>*cX=const_cast<complex<double> *>(X);
    if(X!=Y){
        if(Beta!=0.){
            if(Beta!=1.)scale(Beta,Y,SizY);
            if(Alfa==1.){
                Map<VectorXcd>(Y,SizY).noalias() +=
                        Map<MatrixXcd>(dat->data(),SizY,SizX)
                        * Map<VectorXcd>(cX,SizX);
            }
            else {
                Map<VectorXcd>(Y,SizY).noalias()   +=Alfa*
                        Map<MatrixXcd>(dat->data(),SizY,SizX)
                        * Map<VectorXcd>(cX,SizX);
                // playground for different attempts
                //                matmul(Alfa,dat->data(),cX,SizX,Y,SizY);
                //the following is slower that Eigen:
                //                for(complex<double>*x=cX,*m=dat->data();x!=X+SizX;x++){
                //                    complex<double> aX=*x*Alfa;
                //                    for(complex<double>*y=Y;y!=Y+SizY;y++,m++)
                //                        *y+=*m*aX;
                //                }
            }
        } else {
            if(Alfa==1.)
                Map<VectorXcd>(Y,SizY).noalias() =
                        Map<MatrixXcd>(dat->data(),SizY,SizX)
                        * Map<VectorXcd>(cX,SizX);
            else
                Map<VectorXcd>(Y,SizY).noalias()  =Alfa*
                        Map<MatrixXcd>(dat->data(),SizY,SizX)
                        * Map<VectorXcd>(cX,SizX);
        }

    } else {
        // with aliasing
        if(Beta!=0.){
            if(Beta!=1.)scale(Beta,Y,SizY);
            if(Alfa==1.)
                Map<VectorXcd>(Y,SizY)+=
                        Map<MatrixXcd>(dat->data(),SizY,SizX)
                        * Map<VectorXcd>(cX,SizX);

            else
                Map<VectorXcd>(Y,SizY)+=Alfa*
                        Map<MatrixXcd>(dat->data(),SizY,SizX)
                        * Map<VectorXcd>(cX,SizX);
        } else {
            if(Alfa==1.)
                Map<VectorXcd>(Y,SizY)=
                        Map<MatrixXcd>(dat->data(),SizY,SizX)
                        * Map<VectorXcd>(cX,SizX);

            else
                Map<VectorXcd>(Y,SizY)=Alfa *
                        Map<MatrixXcd>(dat->data(),SizY,SizX)
                        * Map<VectorXcd>(cX,SizX);
        }

    }
}
// multiply from the right
void OperatorZG::axpyTranspose(const std::complex<double>&Alfa, const std::complex<double>* X, unsigned int SizX,
                               const std::complex<double>&Beta, std::complex<double>* Y, unsigned int SizY) const{

    // after some consulting of Eigen advice, Y +=Alfa * Matrix * X seems to be the right thing to do
    if(Alfa==0.){
        scale(Beta,Y,SizY);
        return;
    }

    complex<double>*cX=const_cast<complex<double> *>(X);
    if(X!=Y){
        // NO aliasing
        if(Beta!=0.){
            if(Beta!=1.)scale(Beta,Y,SizY);
            if(Alfa==1.){
                Map<VectorXcd>(Y,SizY).noalias() +=
                        Map<MatrixXcd>(dat->data(),SizX,SizY).transpose()
                        * Map<VectorXcd>(cX,SizX);
            }
            else {
                Map<VectorXcd>(Y,SizY).noalias()   +=Alfa*
                        Map<MatrixXcd>(dat->data(),SizX,SizY).transpose()
                        * Map<VectorXcd>(cX,SizX);
            }
        } else {
            if(Alfa==1.)
                Map<VectorXcd>(Y,SizY).noalias() =
                        Map<MatrixXcd>(dat->data(),SizX,SizY).transpose()
                        * Map<VectorXcd>(cX,SizX);
            else
                Map<VectorXcd>(Y,SizY).noalias()  =Alfa*
                        Map<MatrixXcd>(dat->data(),SizX,SizY).transpose()
                        * Map<VectorXcd>(cX,SizX);
        }

    } else {
        // with aliasing
        if(Beta!=0.){
            if(Beta!=1.)scale(Beta,Y,SizY);
            if(Alfa==1.)
                Map<VectorXcd>(Y,SizY)+=
                        Map<MatrixXcd>(dat->data(),SizX,SizY).transpose()
                        * Map<VectorXcd>(cX,SizX);

            else
                Map<VectorXcd>(Y,SizY)+=Alfa*
                        Map<MatrixXcd>(dat->data(),SizX,SizY).transpose()
                        * Map<VectorXcd>(cX,SizX);
        } else {
            if(Alfa==1.)
                Map<VectorXcd>(Y,SizY)=
                        Map<MatrixXcd>(dat->data(),SizX,SizY).transpose()
                        * Map<VectorXcd>(cX,SizX);

            else
                Map<VectorXcd>(Y,SizY)=Alfa *
                        Map<MatrixXcd>(dat->data(),SizX,SizY).transpose()
                        * Map<VectorXcd>(cX,SizX);
        }

    }
}

void OperatorZG::construct(const Eigen::MatrixXcd & Mat, std::string Kind) {
    
    string hash=Kind+tools::str(Mat.size());

    // CHECK row wise or column wise.///
    vector<complex<double> >mDat(Mat.size());
    unsigned int k=0;
    for(unsigned int n=0;n<Mat.cols();n++)
        for(unsigned int m=0;m<Mat.rows();m++,k++)
            mDat[k]=Mat(m,n);

    dat=addComplex(hash,mDat);
    _rows=Mat.rows();
    _cols=Mat.cols();
    oNorm=Mat.lpNorm<Eigen::Infinity>();
}

void OperatorZG::matmul(const complex<double>& Alfa, complex<double>*M,complex<double>*X,unsigned int SizX, complex<double>*Y, unsigned int SizY) {
    {
        Map<VectorXcd>(Y,SizY).noalias()+=Alfa*Map<MatrixXcd>(M,SizY,SizX)*Map<VectorXcd>(X,SizX);

    }
}
OperatorZG::OperatorZG(const std::vector<std::complex<double> > &mDat, unsigned int Rows, unsigned int Cols, string Kind)
    :OperatorFloor(Rows,Cols,"ZG")
{
    string hash=Kind+tools::str(Rows)+tools::str(Cols);
    dat=addComplex(hash,mDat);
    // the next two lines are not needed?
    _rows=Cols;
    _cols=Rows;
    oNorm=0.;
    for(unsigned int k=0;k<mDat.size();k++)oNorm=max(oNorm,abs(mDat[k]));
}

void OperatorZG::uniqueData(){
    string hash="unique"+tools::str(this);
    dat=addComplex(hash,*dat);
}

long OperatorZG::applyCount() const{
    return _rows*_cols;
}
