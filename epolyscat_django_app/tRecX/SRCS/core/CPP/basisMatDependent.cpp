// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisMatDependent.h"

using namespace std;

#include "basisGrid.h"

#include "algebra.h"
#include "algebraMulti.h"
#include "basisDvr.h"
#include "basisMat1D.h"

/// Mat = <IBas0|<IBas1|...<IBasN| Pot |JBasN>...|JBas1>|JBas0>
/// <br> Pot=factor_string<{}><{}>...<algebraic expression of coordinate names> (we might  this with <{algebraic expression}>)
/// <br>  factor_string must be convertible to Algebra
/// <br> example:   Op=0.5<{}><Function(Alg(X1),Alg(X2))> for the HO potential in cartesian coordinates
/// <br> (see tutorial/110Pot2d.inp,111Pot2dCO2.inp,112ScatterNonSpherical.inp)
void BasisMatDependent::_construct(std::string Op, const Index * IIndex, const Index* JIndex){
    string modOp=modify(Op,IIndex,JIndex);
    if(modOp==Op)return;

    string topAx=modOp.substr(Op.rfind(":")+1,Op.rfind(">")-Op.rfind(":")-1);
    vector<double> grid,weig,jgrid;
    vector<complex<double> >ivals,jvals;
    valsOnQuadgrid(topAx,3,IIndex, grid,weig,ivals);
    valsOnQuadgrid(topAx,3,JIndex,jgrid,weig,jvals);
    for(size_t k=0;k<grid.size();k++)
        if(grid[k]!=jgrid[k])DEVABORT("quadrature grids are inconsistent");

    const BasisDVR * b=dynamic_cast<const BasisDVR*>(IIndex->basis());
    if(b==0)ABORT("only for BasisDVR floors, is: "+IIndex->str());
    // get weights and support points for basis functions
    BasisMat1D one("<1>",IIndex->basis(),JIndex->basis());
    BasisMat1D qqq("<Q>",IIndex->basis(),JIndex->basis());
    if(one.isEmpty() or qqq.isEmpty())ABORT(Str("fail")+one.isEmpty()+qqq.isEmpty());

    vector<string>coor=tools::splitString(tools::stringInBetween(modOp,":",">"),':');
    if(coor.size()!=2)DEVABORT("for now, only 2d, is: "+modOp);

    modOp=tools::stringInBetween(modOp,"<",">");
    vector<complex<double> > pot,qDvr,wDvr;
    for(size_t k=0;k<b->size();k++){
        wDvr.push_back(one.mat()(k,k));
        qDvr.push_back(qqq.mat()(k,k)/wDvr.back());
        pot.push_back(0.);
    }
    shared_ptr<const AlgebraMulti>opAlg(AlgebraMulti::factory(modOp));
    if(opAlg){
        // 2-argument function
        for(size_t l=0;l<grid.size();l++){
            for(size_t k=0;k<pot.size();k++)
                pot[k]+=weig[l]*wDvr[k]*std::conj(ivals[l])*jvals[l]*opAlg->valMulti({grid[l],qDvr[k]});
        }
    }
    else {
        // construct one algebra for each DVR point
        if(Op.find("Q")!=std::string::npos){
            ABORT("cannot use \"Q\" in multi-argument function "+Op
                  +"\n use axis names as variables, here: "+IIndex->coordinates());
        }
        vector<std::shared_ptr<Algebra> > opAlgK;
        //  higher axis name with "Q" (it is assumed that this in most cases will be more different numbers)
        string opQ(Op);
        size_t pos;
        while(string::npos!=(pos=opQ.find(topAx)))opQ.replace(pos,topAx.length(),"Q");
        opQ.replace(opQ.find(":Q>"),3,">");
        opQ=tools::stringInBetween(opQ,"<",">");
        for(size_t k=0;k<qDvr.size();k++){
            // replace present (DVR) axis name with with value
            string op(opQ);
            while(string::npos!=(pos=op.find(IIndex->axisName())))
                op.replace(pos,IIndex->axisName().length(),"("+tools::str(qDvr[k].real(),15)+"+i*("+tools::str(qDvr[k].imag(),15)+"))");
            opAlgK.push_back(std::shared_ptr<Algebra>(new Algebra(op)));
            if(not opAlgK.back()->isAlgebra())ABORT("cannot interprete expression: "+op+"\nInfo: "+Algebra::failures);
        }
        for(size_t l=0;l<grid.size();l++)
            for(size_t k=0;k<pot.size();k++){
                pot[k]+=weig[l]*wDvr[k]*std::conj(ivals[l])*jvals[l]*opAlgK[k]->val(grid[l]);
            }

    }
    _mat=Eigen::Map<Eigen::VectorXcd>(pot.data(),pot.size()).asDiagonal();
}

// append currend axis name to last factor: <{}><term> ->  <{}><term:axNam> (if DVR, can force by <{Int}><term>
std::string BasisMatDependent::modify(const string Op, const Index *IIndex, const Index *JIndex){


    if(Op.find("<{}>")!=Op.find("<") and
            Op.find("<{Int}>")!=Op.find("<") and
            Op.rfind(":")>=Op.find(">"))
        return Op; // not a dependent operator

    if(Op.find(IIndex->axisName(),Op.rfind("<"))==string::npos){
        PrintOutput::warning("no dependence of "+Op+" on present axis "+IIndex->axisName());
        return Op; // no real dependencies
    }

    if(Op.find("<{}>")!=string::npos and IIndex->basis()->isDVR() and  JIndex->basis()->isDVR())return Op; // DVR integrate if possible

    if(IIndex->axisName()!=JIndex->axisName())DEVABORT("not the same axis\n"+IIndex->strNode()+"\n"+JIndex->strNode());

    string s(Op);
    s.insert(s.rfind(">"),":"+IIndex->axisName());
    return s;
}

void BasisMatDependent::valsOnQuadgrid(const std::string TopAx, int Inflate, const Index* Idx,
                                       vector<double> &Grid, vector<double> &Weig, vector<std::complex<double> > &Vals){
    // get dependent function at topAx on "good" quadrature grid
    const Index *iBranch=Idx;
    while(iBranch!=0 and iBranch->axisName()!=TopAx)iBranch=iBranch->parent();
    if(iBranch==0)ABORT("dependence axis "+TopAx+" not found in hierachy: "+Idx->root()->hierarchy());

    const BasisIntegrable *b=iBranch->basis()->integrable();
    //    const BasisDVR *b=dynamic_cast<const BasisDVR*>(iBranch->basis());
    if(b==0)ABORT("not an integrable basis: "+iBranch->basis()->str());
    UseMatrix grid,weig,vals;
    b->quadRule(b->order()*Inflate,grid,weig);
    //    b->dvrRule(grid,weig);
    vals=b->val(grid);


    Grid.clear();
    Weig.clear();
    Vals.clear();
    int nfunc=Idx->index()[iBranch->depth()];
    for(size_t k=0;k<grid.size();k++){
        Grid.push_back(grid(k).real());
        Weig.push_back(weig(k).real());
        Vals.push_back(vals(k,nfunc).complex());
    }

}
