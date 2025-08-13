// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "eigenTools.h"

#include "useMatrix.h"

#include <stdio.h>      /* fopen, fputs, fclose, stderr */
#include <iostream>
#include <fstream>
#include <sstream>
#include "abort.h"
#include "tools.h"
#include "str.h"

namespace EigenTools {

bool isNaN(const Eigen::MatrixXcd & M){
    for(const std::complex<double> *a=M.data();a-M.data()<M.size();a++)
        if(std::isnan(a->real()) or std::isnan(a->imag()))return true;
    return false;
}

bool isIdentity(const Eigen::MatrixXcd & Mat, double Eps){
    return UseMatrix::UseMap(const_cast<Eigen::MatrixXcd*>(&Mat)->data(),Mat.rows(),Mat.cols()).isIdentity(Eps);
}

Eigen::MatrixXcd purge(const Eigen::MatrixXcd &Mat, double Eps){
    Eigen::MatrixXcd m(Mat);
    return purgeInPlace(m,Eps);
}
Eigen::MatrixXcd& purgeInPlace(Eigen::MatrixXcd &Mat, double Eps){
    for(int k=0;k<Mat.size();k++){
        double r=Mat.data()[k].real();
        double i=Mat.data()[k].imag();
        Mat.data()[k]=std::complex<double>(std::abs(r)<Eps?0:r,std::abs(i)<Eps?0:i);
    }
    return Mat;
}

bool isZero(const Eigen::MatrixXcd & Mat, double Eps){
    return UseMatrix::UseMap(const_cast<Eigen::MatrixXcd*>(&Mat)->data(),Mat.rows(),Mat.cols()).isZero(Eps);
}

static bool symmetry(const Eigen::MatrixXcd & Mat,double Eps, int Kind){
    if(Mat.rows()!=Mat.cols())return false;
    for(int k=0;k<Mat.rows();k++){
        double epsSq=Mat.col(k).squaredNorm()*Eps*Eps;
        if(Kind==0){
            if((Mat.row(k)-Mat.col(k).adjoint()).squaredNorm()>epsSq)return false;
        }
        else if(Kind==1){
            if((Mat.row(k)-Mat.col(k).transpose()).squaredNorm()>epsSq)return false;
        }
        else if(Kind==2){
            if((Mat.row(k)+Mat.col(k).transpose()).squaredNorm()>epsSq)return false;
        }
        else
            ABORT(Sstr+"symmetry Kind: 0..s.a., 1..symmetric, 2..antisymmetric, found: "+Kind);
    }
    return true;

}
bool isSelfadjoint(const Eigen::MatrixXcd & Mat, double Eps){return symmetry(Mat,Eps,0);}
bool isSymmetric(const Eigen::MatrixXcd & Mat, double Eps){return symmetry(Mat,Eps,1);}
bool isAntisymmetric(const Eigen::MatrixXcd & Mat, double Eps){return symmetry(Mat,Eps,2);}


std::string str(const Eigen::MatrixXcd &Mat, int Digits, char Print, int Width){

    if(Print=='A'){
        // automatically decide what to print (default)
        double rNorm=Mat.real().lpNorm<Eigen::Infinity>(),iNorm=Mat.imag().lpNorm<Eigen::Infinity>();
        Print='C';
        if(iNorm<rNorm*1e-10)Print='R';
        if(rNorm<iNorm*1e-10)Print='I';
    }

    std::ostringstream oss;
    // column numbering
    int width=Width?Width:Digits+7;

    if(Digits>0)
        for (unsigned int n=0;n<Mat.cols();n++){
            if(n==0){
                oss<<std::setw(width+3);
                if     (Print=='R')oss<<"(real) 0";
                else if(Print=='I')oss<<"(imag) 0";
                else oss<<0;
            }
            else oss<<std::setw(width)<<n;
        }
    else {
        oss<<"   ";
        for (int n=10;n<Mat.cols();n+=10)oss<<"|"<<std::setw(9)<<n;
        oss<<"|";
    }
    oss<<std::endl;

    for (unsigned int m=0;m<Mat.rows();m++){
        // row numbering
        oss<<std::setw(4)<<m;
        // real parts
        for (unsigned int n=0;n<Mat.cols();n++){
            if(Digits>0)oss<<std::setprecision(Digits)<<std::setw(width)
                          <<(Print=='I'?Mat(m,n).imag():Mat(m,n).real());
            else if(Print=='R')oss<<tools::zero(Mat(m,n).real());
            else if(Print=='I')oss<<tools::zero(Mat(m,n).imag());
            else oss<<tools::zero(Mat(m,n));
        }
        oss<<std::endl;

        // imaginary parts
        if(Digits>0 and Print=='C'){
            oss<<std::setw(4)<<" ";
            if(Digits>0)
                for (unsigned int n=0;n<Mat.cols();n++)
                    oss<<std::setprecision(Digits)<<std::setw(width)<<Mat(m,n).imag();
            oss<<std::endl;
        }
    }
    if(Print=='R')oss<<" (real) ";
    else if(Print=='I')oss<<" (imag) ";
    else oss<<" (complex) ";

    oss<<Mat.rows()<<" x "<<Mat.cols()<<std::endl;
    return oss.str();
}

void saveAscii(std::string FileName, const UseMatrix &Mat){
    saveAscii(FileName,Eigen::Map<Eigen::MatrixXcd>(Mat.data(),Mat.rows(),Mat.cols()));
}
void saveAscii(std::string FileName, const Eigen::MatrixXcd &Mat){
    std::ofstream stream(FileName.c_str(),std::ios::out);
    if(not stream.is_open())ABORT("failed to open write file "+FileName);
    stream<<"Info: cd"<<std::endl;
    stream<<Mat.rows()<<" "<<Mat.cols()<<std::endl;
    stream<<std::setprecision(16);
    for(int j=0;j<Mat.cols();j++)
        for(int i=0;i<Mat.rows();i++)
            stream<<i<<" "<<j<<" "<<Mat(i,j)<<std::endl;
    std::cout<<"saveAscii matrix to "+FileName<<std::endl;
    stream.close();
}

void readAscii(std::string FileName, Eigen::MatrixXcd &Mat){
    if(FileName==""){
        std::cout<<"no FileName given, enter here:"<<std::endl;
        std::cin>>FileName;
    }
    if(FileName=="")ABORT("no file name");
    std::ifstream stream(FileName.c_str(),std::ios::in);
    if(not stream.is_open())ABORT("failed to open read file "+FileName);
    std::string line;
    std::getline(stream,line);
    std::getline(stream,line);
    int rows,cols;
    std::stringstream(line)>>rows>>cols;
    Mat=Eigen::MatrixXcd::Zero(rows,cols);
    int i,j;
    while(std::getline(stream,line)){
        std::stringstream(line)>>i>>j;
        std::stringstream(line)>>i>>j>>Mat(i,j);
    }
}

bool compareMatrices(std::string FileA, std::string FileB, double Epsilon){
    Eigen::MatrixXcd A,B,C;
    readAscii(FileA,A);
    readAscii(FileB,B);
    C=A-B;
    if(A.rows()!=B.rows() or A.cols()!=B.cols()){
        std::cout<<"Fail: Matrix dimensions do not match "<<A.rows()<<"X"<<A.cols()<<" vs. "<<B.rows()<<"X"<<B.cols()<<std::endl;
        return false;
    }
    else if(not C.isZero(Epsilon)){
        for(int j=0;j<A.cols();j++)
            for(int i=0;i<A.rows();i++)
                if(std::abs(C(i,j))>std::max(1.,abs(A(i,j))+abs(B(i,j)))*Epsilon)
                    std::cout<<i<<j<<": A="<<A(i,j)<<" B="<<B(i,j)<<std::endl;
        std::cout<<str(C,0)<<std::endl; // show structure of error
        std::cout<<"FAIL: "<<A.rows()<<"x"<<A.cols()<<" matrices from "<<FileA<<" and "<<FileB<<" disagree"<<std::endl;
        return false;
    }

    std::cout<<"OK: "<<A.rows()<<"x"<<A.cols()<<" matrices from "<<FileA<<" and "<<FileB<<" agree"<<std::endl;
    return true;
}

}
