// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisMat1D.h"

#include "basisAbstract.h"
#include "algebra.h"

#include <vector>
using namespace std;

#include "basisIntegrable.h"
#include "basisDvr.h"
#include "basisSub.h"
#include "basisGridQuad.h"
#include "tools.h"
static const double ErrorDouble=DBL_MAX;

/// Mat = <IBas|Op|JBas>, IBas and JBas are the bases at IIndex and JIndex
/// <br> Op=factor_string<operator_string>
/// <br> factor_string must be a valid string for class Algebra
/// <br> operator_string any of <alg>,<d_alg>,<alg_d> and <d_alg_d>
/// with alg a valid string for class Algebra
BasisMat1D::BasisMat1D(std::string Op, const Index * IIndex, const Index* JIndex){
    if(not IIndex->isBottom())return;

    BasisMat1D bm(Op,IIndex->basis(),JIndex->basis());
    _mat=bm._mat;
}

BasisMat1D::BasisMat1D(std::string Op, int Idim, int Jdim){
    if(Op.find("<Id>")==1){
        if(Idim!=Jdim)return;
        _mat=Eigen::MatrixXcd::Identity(Idim,Jdim);
    }
    else if (Op.find(",")<Op.find_first_of(">[({")){
        int i=std::stoi(Op.substr(Op.find('<')+1,Op.find(',')));
        if(i<0 or i>=Idim)return;
        int j=std::stoi(Op.substr(Op.find(',')+1,Op.find('>')));
        if(j<0 or j>=Jdim)return;
        _mat=Eigen::MatrixXcd::Zero(Idim,Jdim);
        _mat(i,j)=1;
    }
    else if (Op.find("<delta[")==0){
        _mat=Eigen::MatrixXcd::Zero(Idim,Jdim);
        int off=int(Algebra::realConstant(tools::stringInBetween(Op,"[","]")));
        if(off>0)for(size_t k=   0;  int(k)<min(Idim,Jdim-off);k++)_mat(k,k+off)=1.;
        else     for(size_t k=-off;  int(k)<min(Idim,Jdim    );k++)_mat(k,k+off)=1.;
    }
    else if (Op.find("<diagonal[")==0){
        _mat=Eigen::MatrixXcd::Zero(Idim,Jdim);
        Algebra a(tools::stringInBetween(Op,"[","]"));
        if(not a.isAlgebra())ABORT("not an algebra string in "+Op);
        for(size_t k=0;int(k)<min(Idim,Jdim);k++)_mat(k,k)=a.val(k);
    }
}

static bool Basis_hasOverlap(const BasisIntegrable* A,const BasisIntegrable*B){
    return B->lowBound()<A->upBound() and A->lowBound()<B->upBound();
}

BasisMat1D::BasisMat1D(string Op, const BasisAbstract *IBas, const BasisAbstract *JBas){

    //HACK until better solution======================================================

    if(Op=="<Id>" and *IBas==*JBas){
        _mat=Eigen::MatrixXcd::Identity(IBas->size(),JBas->size());
        return;
    }

    const BasisAbstract * iBas=BasisSub::superBas(IBas);
    const BasisAbstract * jBas=BasisSub::superBas(JBas);

    if(Op.find("<GridWeight")==0){
        std::unique_ptr<Algebra> alg;
        if(not(*iBas==*jBas))DEVABORT("GridWeight "+jBas->str()+" != "+JBas->str());
        if(dynamic_cast<const BasisGridQuad*>(iBas)){
            const BasisGridQuad* g=dynamic_cast<const BasisGridQuad*>(iBas);
            _mat=Eigen::MatrixXcd::Zero(g->size(),g->size());
            for(size_t k=0;k<g->size();k++)_mat(k,k)=g->weights()[k];
        }
        else if(dynamic_cast<const BasisGrid*>(iBas)){
            const BasisGrid* g=dynamic_cast<const BasisGrid*>(iBas);
            _mat=Eigen::MatrixXcd::Zero(g->size(),g->size());
            if(g->size()==1)
                _mat(0,0)=1.;
            else if(g->size()==2){
                _mat(0,0)=0.5*(g->mesh()[1]-g->mesh()[0]);
                _mat(1,1)=_mat(0,0);
            }
            else {
                for(size_t k=1;k<g->size()-1;k++)_mat(k,k)=0.5*(g->mesh()[k+1]-g->mesh()[k-1]);
                _mat(0,0)=_mat(1,1);
                _mat(g->size()-1,g->size()-1)=_mat(g->size()-2,g->size()-2);
            }
        }
        else
            DEVABORT("cannot do "+Op+" for basis"+IBas->str());
        if(Op.find("<GridWeight*")==0){
           Algebra alg(tools::stringInBetween(Op,"<GridWeight*",">",true));
           if(not alg.isAlgebra())ABORT("illegal algebra in "+Op);
           std::vector<double> mesh(dynamic_cast<const BasisGrid*>(iBas)->mesh());
           for(size_t k=0;k<mesh.size();k++)_mat(k,k)*=alg.val(mesh[k]);
        }
        _mat=BasisSub::subMatrix(_mat,BasisSub::subset(IBas),BasisSub::subset(JBas));
        return;
    }


    if(iBas->isGrid() or jBas->isGrid()
            or iBas->isIndex() or jBas->isIndex() or iBas->name()=="CIion" or jBas->name()=="CIion"){
        if(Op=="<1>"){
            _mat=BasisSub::subMatrix(Eigen::MatrixXcd::Identity(iBas->size(),jBas->size()),
                                     BasisSub::subset(IBas),BasisSub::subset(JBas));
            return;
        }
        if(tools::findFirstOutsideBrackets(Op,",","[","]")!=string::npos){
            // operators of the form <i,j> (for grid- or index-type basis)
            vector<string> ij=tools::splitString(Op,',');
            if(ij.size()==2){
                _mat=Eigen::MatrixXcd::Zero(iBas->size(),iBas->size());
                _mat(std::stoi(ij[0].substr(1)),std::stoi(ij[1].substr(0,ij[1].find('>'))))=1.;
                return;
            }
        }
    }
    if(iBas->isIndex() and jBas->isIndex()){
        BasisMat1D m(Op,iBas->size(),jBas->size());
        if(not m.isEmpty()){
            _mat=BasisSub::subMatrix(m.mat(),BasisSub::subset(IBas),BasisSub::subset(JBas));
            return;
        }
    }
    //==================================================================================
    // integrable bases and local operators only
    const BasisIntegrable * bi,*bj;
    bi=dynamic_cast<const BasisIntegrable*>(iBas);
    bj=dynamic_cast<const BasisIntegrable*>(jBas);
    if(bi==0 or bj==0)return;
    if(not Basis_hasOverlap(bi,bj))return;


    // get exact or dvr quadrature
    UseMatrix quadX,quadW;

    if(not (*bi==*bj)){
        // bases differ
        //HACK - overkill gauss-legendre quadrature rule
        double a=std::max(bi->lowBound(),bj->lowBound());
        double b=std::min(bi->upBound(),bj->upBound());
        OrthogonalLegendre leg;
        std::vector<double>q,w;
        leg.quadratureGauss(std::max(bi->order(),bj->order())+10,q,w);
        quadX=UseMatrix(q.size(),1);
        quadW=UseMatrix(w.size(),1);
        for(size_t k=0;k<q.size();k++){
            quadW(k)= 0.5*(b-a)* w[k];
            quadX(k)=(0.5*(b-a)*(q[k]+1.))+a;
        }
    }
    else
        if(iBas->isDVR() and jBas->isDVR()){
        const BasisDVR* dv;
        if(0!=(bi=dv=dynamic_cast<const BasisDVR*>(iBas))){
            dv->dvrRule(quadX,quadW);
        }
        else DEVABORT("isDVR, but not BasisDVR: "+iBas->str());
    }
    else if (0!=(bi=dynamic_cast<const BasisIntegrable*>(iBas)) and 0!=(bj=dynamic_cast<const BasisIntegrable*>(jBas)))
        bi->quadRule(max(bi->order()+3,bj->order()+3),quadX,quadW);
    else
        return;


    bj=dynamic_cast<const BasisIntegrable*>(jBas);
    bi->jacobian()->operator ()(quadX,quadW);

    // functions may have discontinuities or  singularities at the interval boundaries
    // --- evaluate slightly away from boundaries ---
    // epsilon for moving away from boundary
    double epsLo=1.e-14,epsUp=1.e-14;
    if(bi->lowBound()>=-BasisIntegrable::infty)epsLo=epsLo+1.e-14*abs(bi->lowBound());
    if(bi->upBound() <= BasisIntegrable::infty)epsUp=epsUp+1.e-14*abs(bi->upBound());

    // this should go to basis integrable, maybe into the complexScaling() function
    // make sure complex scaling radius is at boundary
    double eps=min(1.e-8,(bi->upBound()-bi->lowBound())*1.e-12);
    if(bi->complexScaling()->r0up()>bi->lowBound()+eps and bi->complexScaling()->r0up()<bi->upBound()-eps)
        ABORT(Str("complex scaling radius must be on element boundary, is inside: lb=")+
              bi->lowBound()+"< r0="+bi->complexScaling()->r0up()+"< ub="+bi->upBound());

    UseMatrix funcArg;
    funcArg=quadX;
    for(size_t k=0;k<quadX.size();k++){
        funcArg(k)=bi->complexScaling()->xScaled(quadX(k).real());
        quadW(k)*=bi->complexScaling()->etaX(0.5*(bi->lowBound()+bi->upBound()));
        if(abs(funcArg(k).real()-bi->lowBound())<=epsLo)funcArg(k)=complex<double>(bi->lowBound()+epsLo,funcArg(k).imag());
        if(abs(funcArg(k).real()-bi->upBound()) <=epsUp)funcArg(k)=complex<double>(bi->upBound() -epsUp,funcArg(k).imag());
    }

    complex<double> preFac;
    string fac=Op.substr(0,Op.find("<"));
    if(fac=="")fac="1";
    Algebra facAlg(fac);
    if(not facAlg.isAlgebra())return;
    preFac=facAlg.val(0.);

    string opFunc=tools::stringInBetween(Op,"<",">");
    if(opFunc.find("d_")==0)opFunc=opFunc.substr(2);
    if(opFunc.find("_d")==opFunc.length()-2)opFunc=opFunc.substr(0,opFunc.length()-2);

    Algebra funcAlg(opFunc);
    if(not funcAlg.isAlgebra())return; // not a legal algebra string

    bool warnNonAnalytic=false;
    for(auto sing: funcAlg.nonAnalyticQ())
        warnNonAnalytic=warnNonAnalytic or (bi->lowBound()+epsLo<sing.real() and (sing.real()<bi->upBound()-epsUp));
    if(warnNonAnalytic){
        PrintOutput::warning(Sstr+funcAlg.definition()+" has non-analytic point in quadrature interval ["
                             +bi->lowBound()+","+bi->upBound()+"]",1,0,
                             Sstr+"Exact numerical integration needs integrands that are analytic within every finite element."
                                  "\n Matrix elements will be discretization dependent and non-systematic errors may appear"+
                             +"\n    - errors may be acceptably low, but convergence with order may be inconsistent"
                             +"\n    - for consistent convergence place element boundaries on non-analytic points");
    }

    // we should check compatibility
    Eigen::MatrixXcd iVals(valDer(quadX,bi,Op.find("<d_")!=string::npos));
    Eigen::MatrixXcd jVals(valDer(quadX,bj,Op.find("_d>")!=string::npos));

    Eigen::VectorXcd funcWeig=Eigen::VectorXcd::Zero(quadX.size());
    for(size_t k=0;k<quadX.size();k++){
        if(iVals.row(k).isZero() or jVals.row(k).isZero())continue; // skip zero entries (mostly for avoiding singulariets at r=0)
        funcWeig(k)=funcAlg.val(funcArg(k).complex())*quadW(k).complex();
    }

    _mat=iVals.adjoint()*funcWeig.asDiagonal()*jVals*preFac;
    _mat.purge();

    // complex scaling
    _mat/=std::pow(bi->complexScaling()->etaX(0.5*(bi->lowBound()+bi->upBound())),tools::subStringCount(Op,"_d>")+tools::subStringCount(Op,"<d_"));
    _mat=BasisSub::subMatrix(_mat,BasisSub::subset(IBas),BasisSub::subset(JBas));

}

Eigen::MatrixXcd BasisMat1D::valDer(const UseMatrix &X, const BasisIntegrable *Bas, int Derivative) const
{
    vector<complex<double> > x,val,dum;
    for(size_t k=0;k<X.size();k++)x.push_back(X.data()[k]);
    if(     Derivative==0)Bas->valDer(x,val,dum);
    else if(Derivative==1)Bas->valDer(x,dum,val);
    vector<int> allIdx;
    for(size_t k=0;k<X.size();k++)allIdx.push_back(k);

    return BasisSub::subMatrix(Eigen::Map<Eigen::MatrixXcd>(val.data(),x.size(),val.size()/x.size()),
                               allIdx,BasisSub::subset(Bas));
}
