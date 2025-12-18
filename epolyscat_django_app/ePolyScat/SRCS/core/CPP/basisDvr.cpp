// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisDvr.h"

#include "readInput.h"
#include "printOutput.h"
#include "useMatrix.h"
#include "basisSetDef.h"
#include "tRecXchecks.h"
#include "basisGridQuad.h"

using namespace std;
static bool abs_compare(std::complex<double> a, std::complex<double> b){return (std::abs(a) < std::abs(b));}
bool BasisDVR::femDVR=true;

bool BasisDVR::isDVR() const {
    return _isDVR and
            not (lowBound()<=0. and 0.<=upBound()
                 and  ReadInput::main.flag("DEBUGfem","run FEM basis - compute exact integrals everywhere"));
}

static void getBasicDvr(std::shared_ptr<OrthogonalPolynomial>Opol,
                        const BasisSetDef & Def,std::vector<double> & QuadQ,std::vector<double> &QuadP){

    bool enforceZeroLow=Def.coor.zeroLow and Def.lowBound()>-BasisDVR::infty; // lower bound = -infty implies zero, do not enforce
    bool enforceZeroUp =Def.coor.zeroUp  and Def.upBound()<  BasisDVR::infty; // upper bound =  infty implies zero, do not enforce

    bool leftPoint=not Def.first or enforceZeroLow;
    bool rightPoint=not Def.last or enforceZeroUp;

    if(Def.scale<0.)swap(leftPoint,rightPoint); // interval is reverted

    if(leftPoint and rightPoint)
        Opol->quadratureLobatto(Def.order,QuadQ,QuadP);
    else if(leftPoint)
        Opol->quadratureRadauLeft(Def.order,QuadQ,QuadP);
    else if(rightPoint)
        Opol->quadratureRadauRight(Def.order,QuadQ,QuadP);
    else
        Opol->quadratureGauss(Def.order,QuadQ,QuadP); // no need for quadrature points at ends
}

BasisDVR::BasisDVR(const BasisSetDef &Def):BasisIntegrable(Def.lowBound(),Def.upBound()),_isDVR(true){
    if(Def.funcs=="polynomial"){
        _opol=std::shared_ptr<OrthogonalPolynomial>(new OrthogonalLegendre());
        _scale=Def.scale/(_opol->upperBoundary()-_opol->lowerBoundary());
        _shiftQ=_opol->lowerBoundary();
    }
    else if(Def.funcs.find("jacobi[")!=std::string::npos){
        std::vector<double>ab=tools::string_to_doubleVec("{"+tools::stringInBetween(Def.funcs,"[","]")+"}");
        if(ab.size()!=2)ABORT("need exactly jacobi[a,b], got: "+Def.funcs);
        _opol=std::shared_ptr<OrthogonalPolynomial>(new OrthogonalJacobi(ab[0],ab[1]));
        _scale=Def.scale/(_opol->upperBoundary()-_opol->lowerBoundary());
        _shiftQ=_opol->lowerBoundary();
    }
    else if(Def.funcs=="custom") {
        // [Def.shift , Def.shift+Def.scale] in X coords
        // [-1;+1] in Q coords
        // w(X) = X
        // w(Q) = (Q+1)/2*Def.scale + Def.shift
        //      = Def.scale/2 * (Q + 1+Def.shift*2/Def.scale)
        // first factor will be multiplied into the basis later anyway
        _opol=std::shared_ptr<OrthogonalPolynomial>(new OrthogonalCustom(Def.scale/2, -1.-Def.shift*2./Def.scale));
        _scale=Def.scale/(_opol->upperBoundary()-_opol->lowerBoundary());
        _shiftQ=_opol->lowerBoundary();
    }
    else if(Def.funcs.find("polExp")==0){
        _opol=std::shared_ptr<OrthogonalPolynomial>(new OrthogonalLaguerre());
        _scale=Def.scale;
        _shiftQ=_opol->lowerBoundary();
    }
    else  ABORT("no DVR defined for, is: "+Def.funcs);

    _name="DVR: "+_opol->name();
    _shiftX=Def.shift;
    _lowBound=Def.lowBound();
    _upBound=Def.upBound();
    if(abs(_lowBound)<abs(_scale)*1.e-12)_lowBound=0.;
    if(abs(_upBound)<abs(_scale)*1.e-12)_upBound=0.;

    // if there are Dirichlet conditions, we want a node at boundary

    // get Gauss, Lobatto, or Radau quadrature, depending on boundary conditions and position in axis
    getBasicDvr(_opol,Def,_dvrX,_dvrW);

    for(size_t k=0;k<_dvrX.size();k++){
        _dvrW[k]*=std::abs(_scale)/_opol->weight(_dvrX[k]);
        _dvrX[k]=xFromQ(_dvrX[k]);
    }
    if(_scale<0. or (Def.first and not Def.last and _dvrX[0]==_lowBound and _dvrX.back()!=_upBound)){
        // reverted Radau:
        std::reverse(_dvrW.begin(),_dvrW.end());
        std::reverse(_dvrX.begin(),_dvrX.end());
        if(_scale>0.)
            for(size_t k=0;k<_dvrX.size();k++)_dvrX[k]=_lowBound+_upBound-_dvrX[k];
    }

    if(Def.coor.jaco==Coordinate::J_one)
        _jac=std::shared_ptr<Jacobian>(new Jacobian1(0.));
    else if(Def.coor.jaco==Coordinate::J_val)
        _jac=std::shared_ptr<Jacobian>(new JacobianQ(0.));

    _comSca=std::shared_ptr<const ComplexScaling>(new ComplexScaling(Def.comSca));

    // snap quadrature points at boundary to exact boundary value
    snapBoundary(_dvrX[0],    _lowBound,_upBound);
    snapBoundary(_dvrX.back(),_lowBound,_upBound);

    // enforce boundary conditions and check
    _nBeg=0;
    _size=_dvrX.size();
    vector<complex<double> > x,v,d;
    for(size_t k=0;k<_dvrX.size();k++)x.push_back(_dvrX[k]);
    valDer(x,v,d,false);
    double epsV=1.e-12*abs(*std::max_element(v.begin(),v.end(),abs_compare));
    for(size_t n=0,kn=0;n<_dvrX.size();n++,kn+=_dvrX.size()+1){
        if(abs(v[kn])<epsV*sqrt(_opol->weight(qFromX(_dvrX[n])))){
            PrintOutput::warning(Str("basis evaluates to zero at DVR point n=")+n+"x="+_dvrX[n]+"v="+v[kn]+kn+v.size());
        }

        if(Def.first and n==0 and _dvrX[0]==_lowBound and Def.coor.zeroLow){
            // boundary condition=0 at left end of axis
            _nBeg=1;
            _size--;}
        else if(Def.last and n==_dvrX.size()-1 and _dvrX.back()==_upBound and Def.coor.zeroUp){
            // boundary condition=0 at right end of axis
            _size--;
        }
        //        else if(abs(v[kn])<epsV*sqrt(_opol->weight(qFromX(_dvrX[n]))))
        //            ABORT(Str("basis=")+abs(v[kn])+"at DVR point, but not axis boundary and zero boundary condition, x="+_dvrX[n]);
    }
    setDVR(not Def.exactIntegral);
}
std::vector<double> BasisDVR::valNodes() const {
    // actual values at collocation points (the safe way)
    std::vector<std::complex<double>>xc,vc,dc;
    std::vector<double> _valNodes;
    for(double x: _dvrX)xc.push_back(x);
    valDer(xc,vc,dc,false);
    for(size_t k=0;k<size();k++)_valNodes.push_back(vc[_nBeg+k+k*xc.size()].real());
    return _valNodes;
}

void BasisDVR::dvrRule(UseMatrix & QuadX, UseMatrix & QuadW) const{
    QuadX=UseMatrix(_dvrX.size(),1);
    QuadW=UseMatrix(_dvrW.size(),1);
    for(size_t k=0;k<_dvrX.size();k++){
        QuadX.data()[k]=_dvrX[k];
        QuadW.data()[k]=_dvrW[k];
    }
}

void BasisDVR::quadRule(int N,std::vector<double>&QuadX,std::vector<double>&QuadW) const
{
    // this, eventually, should be moved to BasisIntegrable
//    if(_opol->name()!="legendre" && _opol->name().find("laguerre")!= 0 && _opol->name().find("jacobi")!= 0)DEVABORT("needs implementation");
    _opol->quadrature(N,QuadX,QuadW);
    for(size_t k=0;k<QuadX.size();k++){
        QuadW[k]/=_opol->weight(QuadX[k]);
    }
    for(size_t k=0;k<QuadX.size();k++){
        QuadX[k]=xFromQ(QuadX[k]);
        QuadW[k]=QuadW[k]*std::abs(_scale);
    }
}

void BasisDVR::valDer(const std::vector<std::complex<double> >&X,
                      std::vector<std::complex<double> >&Val,std::vector<std::complex<double> >&Der, bool ZeroOutside) const{
    // see also tsurff.pdf

    // starting values sqrt(v(q(x))), d/dx sqrt(v(q(x))
    std::vector<double> val0,der0;
    for(size_t i=0;i<X.size();i++){
        val0.push_back(    sqrt(   _opol->weight(qFromX(X[i].real()))));
        der0.push_back(0.5*sqrt(1./_opol->weight(qFromX(X[i].real())))
                       *_opol->derWeight(qFromX(X[i].real()))/_scale);
    }

    Val.clear();
    Der.clear();
    // column index n...function numbers b_n, row index i...x-values x[i]
    for(int n=_nBeg;n<_nBeg+_size;n++){
        // tabulate divisions
        std::vector<double>qDiffN(_dvrX.size(),0.);
        for(int k=0;k<int(_dvrX.size());k++)
            if(k!=n)qDiffN[k]=1./(_dvrX[n]-_dvrX[k]);

        // get norm of n'th function; on margin, force value=1, else ||b_n||=1
        double normN=1./_opol->weight(qFromX(_dvrX[n]));
        if(not (   abs(_lowBound-_dvrX[n])<abs(_scale)*1e-10
                   or abs(_upBound-_dvrX[n])<abs(_scale)*1e-10))normN/=_dvrW[n];
        normN=sqrt(normN);

        for(size_t i=0;i<X.size();i++){
            if(ZeroOutside and (X[i].real()<_lowBound or X[i].real()>_upBound)){
                Val.push_back(0.);
                Der.push_back(0.);
            }
            else {
                Val.push_back(val0[i]*normN);
                Der.push_back(der0[i]*normN);
                for(size_t k=0;k<_dvrX.size();k++){
                    if(int(k)==n)continue;
                    Der.back()=(Der.back()*(X[i].real()-_dvrX[k])+Val.back())*qDiffN[k];
                    Val.back()=(Val.back()*(X[i].real()-_dvrX[k])           )*qDiffN[k];
                }
            }
        }
    }
}

std::string BasisDVR::str(int Level) const{
    std::string res=Str(name(),"")+" ["+_lowBound+","+_upBound+"] ("+_shiftX+","+_scale+") "+_size+"["+_dvrX.size()+"]";
    return isDVR()?res:res+" noDVR";
}

std::string BasisDVR::strDefinition() const {
    //NOTE: much of this info is redundand, but reconstruction from minimal info is complex, needs care
    const BasisGridQuad * grid=BasisGridQuad::factory(_dvrX,_dvrW); // this is redundand, but the easy solution for now
    std::string res=Str(BasisIntegrable::strDefinition(),",")+_shiftX+_scale+_nBeg+("["+grid->strDefinition()+"]");
    return _isDVR?res:res+", noDVR";
}
BasisDVR::BasisDVR(const string &Def):BasisIntegrable(0.,0.),_isDVR(true){
    if(Def.substr(0,4)!="DVR:")DEVABORT("not a DVR definition, is "+Def);

    // split outside the brackets
    std::vector<std::string> def(tools::splitString(Def.substr(4),',',"[(","])"));
    if(def[0]=="legendre"){
        _opol.reset(new OrthogonalLegendre());
        _shiftQ=_opol->lowerBoundary();
    }
    else if(def[0].find("laguerre(")==0){
        _opol.reset(new OrthogonalLaguerre(tools::string_to_int(tools::stringInBetween(def[0],"(",")"))));
        _shiftQ=_opol->lowerBoundary();
    }
    else if(def[0].find("jacobi[")==0){
        std::vector<double> ab=tools::string_to_doubleVec("{"+tools::stringInBetween(def[0],"[","]")+"}");
        _opol.reset(new OrthogonalJacobi(ab[0],ab[1]));
        _shiftQ=_opol->lowerBoundary();
    }
    else{
        DEVABORT("orthogonal polynomial "+def[0]+" not implemented for construction from string\nDef = "+Def);
    }
    // BasisAbstract
    _name="DVR:"+_opol->name();
    _size  =tools::string_to_int(def[1]);

    // BasisIntegrable
    _dvrX.resize(tools::string_to_int(def[2])); // order
    _lowBound=tools::string_to_double(def[3]);
    _upBound=tools::string_to_double(def[4]);
    _comSca.reset(new ComplexScaling(tools::stringInBetween(def[5],"[","]")));
    _jac.reset(Jacobian::factory(tools::cropString(def[6]),_comSca->etaX(_lowBound/2.+_upBound/2.)));

    // DVR-specific
    _shiftX=tools::string_to_double(def[7]);
    _scale =tools::string_to_double(def[8]);
    _nBeg  =tools::string_to_int(def[9]);

        const BasisGridQuad* grid=BasisGridQuad::factory(tools::stringInBetween(def[10],"[","]"));
        _dvrX=grid->mesh();
        _dvrW=grid->weights();
    if(def.size()>11)_isDVR=(def[11]=="noDVR");

    if(abs(_lowBound)<abs(_scale)*1.e-12)_lowBound=0.;
    if(abs(_upBound)<abs(_scale)*1.e-12)_upBound=0.;

}

/// smarter float comparison (as roundoff can be ignored)
bool compareFloat(const std::vector<double> A, const std::vector<double> B, double Eps){
    // fast comparison first;
    if(A.size()!=B.size())return false;
    for(size_t k=0;k<A.size();k++)
        if(abs(A[k]-B[k])>Eps)return false;
    return true;
}


bool BasisDVR::operator==(const BasisAbstract & Other) const {
    if(this==&Other)return true;
    const BasisDVR* other;
    if(0==(other=dynamic_cast<const BasisDVR*>(&Other)))return false;
    if(not compareFloat(_dvrX,other->_dvrX,
                        1.e-12*(max(_dvrX.back(),other->_dvrX.back())-min(_dvrX.front(),other->_dvrX.front()))))
    {
        return false;
    }
    if(not compareFloat(_dvrW,other->_dvrW,1.e-12*(_dvrX.back()-_dvrX.front())))return false;
    return _opol->name()==other->_opol->name();
}
