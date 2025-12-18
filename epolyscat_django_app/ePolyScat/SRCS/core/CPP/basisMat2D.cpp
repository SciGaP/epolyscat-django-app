// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisMat2D.h"

#include "tools.h"
#include "basisAbstract.h"

#include <vector>
using namespace std;

#include "basisMat1D.h"

/// Mat = <IBas0|Op|JBas0><IBas1|Op|JBas1>, IBasN and JBasN are the bases at IIndex and JIndex
/// <br> Op=factor_string<operator_string>
/// <br>  factor_string must be convertible to Algebra
/// <br> operator_string any of <alg>,<d_alg>,<alg_d> and <d_alg_d>
/// with alg convertible to Algebra
BasisMat2D::BasisMat2D(std::string Op, const Index * IIndex, const Index* JIndex){
    if(IIndex->heightAboveBottom()!=1)return;
    if(Op.find("<{}>")!=string::npos)return; // correlated, multiplicative factor

    if(tools::subStringCount(Op,"><")!=1)ABORT("need operator string format <Op1><Op2>, got: "+Op);
    string op1=Op.substr(0,Op.find("><")+1);
    string op2=Op.substr(Op.find("><")+1);
    if(((op1.find("<d_") != string::npos or op1.find("Q") != string::npos) and IIndex->basis()->isGrid())
            or ((op1.find("_d>") != string::npos or  op1.find("Q") != string::npos)and JIndex->basis()->isGrid()) or
            ((op2.find("<d_") != string::npos or op2.find("Q") != string::npos) and IIndex->child(0)->basis()->isGrid()) or
            ((op2.find("_d>") != string::npos or op2.find("Q") != string::npos) and JIndex->child(0)->basis()->isGrid())){
        _mat=Eigen::MatrixXcd::Zero(IIndex->sizeStored(),JIndex->sizeStored());
        return;
    }
    BasisMat1D bm(op1,IIndex->basis(),JIndex->basis());
    if(bm.isEmpty())return;
    _mat0=bm.mat();
    if(not IIndex->subEquivalent() or not JIndex->subEquivalent()){
        // basis is not tensor product
        _mat=Eigen::MatrixXcd::Zero(IIndex->sizeStored(),JIndex->sizeStored());
        for(int j=0;j<_mat0.cols();j++)
            for(int i=0;i<_mat0.rows();i++){
                if(_mat0(i,j)!=0.){
                    BasisMat1D bn(op2,IIndex->child(i)->basis(),JIndex->child(j)->basis());
                    if(bn.isEmpty())return;
                    _mat.block(i*bn.mat().rows(),j*bn.mat().cols(),bn.mat().rows(),bn.mat().cols())=bn.mat()*_mat0(i,j);
                }
            }
    }
    else {
        BasisMat1D bn(op2,IIndex->child(0)->basis(),JIndex->child(0)->basis());
        if(bn.isEmpty())return;
        _mat1=bn.mat();
    }
}

const std::vector<const Eigen::MatrixXcd*> BasisMat2D::mats() const {
    if(_mat0.size()!=0)
        return std::vector<const Eigen::MatrixXcd*>(2)={&_mat0,&_mat1};
    else
        return std::vector<const Eigen::MatrixXcd*>(1,&_mat);
}
