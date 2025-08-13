// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisMatMulti.h"

#include "tools.h"
#include "algebra.h"
#include "basisAbstract.h"
#include "basisDvr.h"

#include <vector>
using namespace std;

#include "basisMat1D.h"

/// Mat = <IBas0|<IBas1|...<IBasN| Pot |JBasN>...|JBas1>|JBas0>
/// <br> Pot=factor_string<{}><{}>...<algebraic expression of coordinate names> (we might replace this with <{algebraic expression}>)
/// <br>  factor_string must be convertible to Algebra
/// <br> example:   Op=0.5<{}><{}><X*X+Y*Y+Z*Z> for the HO potential in cartesian coordinates
void BasisMatMulti::_construct(std::string Op, const Index * IIndex, const Index* JIndex){

     if(Op.find("<{}>")==string::npos)return; // not a multi-dimensional multiplication

    // do the basic checks here...

    string fac=Op.substr(0,Op.find("<"));
    if(fac=="")fac="1";
    Algebra facAlg(fac);
    if(not facAlg.isAlgebra())return;


    // recursively construct values
    std::vector<complex<double> >diag;
    construct(Op.substr(Op.find("<")),IIndex,JIndex,diag);

    Eigen::Map<Eigen::VectorXcd> d(diag.data(),diag.size());
    d*=facAlg.val(0.);
    _mat=d.asDiagonal();
}

void BasisMatMulti::construct(std::string Op, const Index * IIndex, const Index* JIndex, std::vector<std::complex<double> > & Diag){

    // for speed, these checks should be moved to top level
    if(IIndex->axisName()!=JIndex->axisName())
        ABORT(Str("multiplicative operators only between equal coordinates, found:")+IIndex->axisName()+JIndex->axisName()+"for"+Op);
    if(not IIndex->basis()->isDVR())ABORT("only for DVR basis, found "+IIndex->basis()->str());
    if(not (IIndex->basis()->operator==(*JIndex->basis())))
        ABORT(Str("multiplicative operators only for identical basis, found:")+IIndex->basis()->str()+JIndex->basis()->str());
    const BasisDVR* dvr=dynamic_cast<const BasisDVR*>(JIndex->basis());
    if(dvr==0)ABORT("cast to dvr failed");

    if(Op.find("<{}>")==0){
        // complex scaling for weight

        // get weights and support points for basis functions
        BasisMat1D one("<1>",IIndex->basis(),JIndex->basis());
        BasisMat1D qqq("<Q>",IIndex->basis(),JIndex->basis());
        if(one.isEmpty() or qqq.isEmpty())ABORT(Str("fail")+one.isEmpty()+qqq.isEmpty());

        for(size_t k=0;k<IIndex->basis()->size();k++)
        {
            std::complex<double> weigJacValsq=one.mat()(k,k);
            complex<double> q=(qqq.mat()(k,k)/weigJacValsq);

            // substitute coordinate with q value
            string op(Op.substr(Op.find("<{}>")+4));
            size_t pos;
            // Algebra accepts complex constant only as (constR+i*(constI)):
            while(string::npos!=(pos=op.find(IIndex->axisName())))
                op.replace(pos,IIndex->axisName().length(),"("+tools::str(q.real(),15)+"+i*("+tools::str(q.imag(),15)+"))");

            int newDiag=Diag.size(); // begin of new portion
            construct(op,IIndex->child(k),JIndex->child(k),Diag);
            for(size_t k=newDiag;k<Diag.size();k++)Diag[k]*=weigJacValsq;
        }
    }
    else{
        // last factor should hold the algebraic expression
        size_t pos;
        while(string::npos!=(pos=Op.find(IIndex->axisName())))Op.replace(pos,IIndex->axisName().length(),"Q");
        BasisMat1D bm(Op,IIndex->basis(),JIndex->basis());
        if(bm.isEmpty())ABORT("failed to evaluate: "+Op);
        for(int k=0;k<bm.mat().rows();k++)Diag.push_back(bm.mat()(k,k));
    }
}


//const UseMatrix BasisMatMulti::useMat() const {
//    return UseMatrix::UseMap(const_cast<BasisMatMulti*>(this)->_mat.data(),_mat.rows(),_mat.cols());
//}

//const std::vector<const Eigen::MatrixXcd*> BasisMatMulti::mats() const {
//    return std::vector<const Eigen::MatrixXcd*>(1,&_mat);
//}
