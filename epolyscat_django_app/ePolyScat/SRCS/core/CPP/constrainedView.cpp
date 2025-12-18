// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "constrainedView.h"
#include "index.h"
#include "printOutput.h"
#include "coefficients.h"
#include "operatorTree.h"
#include "operatorFloor.h"
#include "algebra.h"
#include "str.h"

#include "parallelOperator.h"
#include "basisIntegrable.h"

using namespace std;

const OperatorAbstract * ConstrainedView::factory(OperatorTree *Op, string Constraint){

    if(Constraint=="")return Op;

    if(Constraint.find("ZeroBoundary=")!=string::npos){
        OperatorTree* op=dynamic_cast<OperatorTree*>(Op);
        if(op==0)ABORT("Constraint "+Constraint+" only for OperatorTree");

        vector<vector<double> > zeros(1);
        zeros[0].push_back(Algebra(Constraint.substr(Constraint.find("ZeroBoundary=")+13)).val(0.).real());
        if(Op->iIndex->hierarchy().find("Rn")==string::npos)ABORT("Constraint::ZeroBoundary only for radial coordinates");
        vector<string>ax(1,"Rn");
        return new ConstrainedZero(op,ax,zeros);
//        return new ConstrainedZeroDerivative(op,ax,zeros);
    }

    ABORT("Constraint not defined: "+Constraint);
    return 0;
}

std::string ConstrainedView::strNode() const {
    Str s("","");
    if(childSize()==0)s=s+"f";
    s=(s+"<"+iIndex->index()+"|"+jIndex->index()+"> ("+iIndex->sizeCompute()+"x"+jIndex->sizeCompute()+") ["+childSize()+"] "+iIndex->axisName());
    if(iIndex->axisName()!=jIndex->axisName())s=s+"-"+jIndex->axisName();
    return std::move(s);
}

// build tree
ConstrainedZero::ConstrainedZero(OperatorTree *Op, std::vector<std::string> Axes, const std::vector<std::vector<double> > &Zeros)
{
    _op=Op;
    name = Op->name + " (zero-constrained)";
    iIndex=Op->iIndex;
    jIndex=Op->jIndex;
    if(Op->isLeaf()){
        vector<unsigned int> leftZeros,rightZeros;
        for(int k=0;k<Axes.size();k++){
            addZeros(Axes[k],Zeros[k],iIndex,leftZeros);
            addZeros(Axes[k],Zeros[k],jIndex,rightZeros);
        }
        if((leftZeros.size()>0 or rightZeros.size()>0))
            oFloor=new Floor(Op->floor(),leftZeros,rightZeros);
        else
        {
            _view=true;
            oFloor=const_cast<OperatorFloor*>(Op->floor());
        }
    }
    for(int k=0;k<Op->childSize();k++)
        childAdd(new ConstrainedZero(Op->child(k),Axes,Zeros));
}

void ConstrainedZero::addZeros(std::string Axis, const std::vector<double> & Zeros, const Index *Idx, std::vector<unsigned int> &IZeros)
{
    if(Idx->heightAboveBottom()>0)ABORT("for now only 1-d constraints");

    if(Idx->axisName()==Axis and Idx->basis()->integrable()!=0){
        double lb=Idx->basis()->integrable()->lowBound();
        double ub=Idx->basis()->integrable()->upBound();
        for(int l=0;l<Zeros.size();l++){
            if(lb+1.e-12<Zeros[l] and Zeros[l]<ub-1.e-12)
                ABORT(Str("Zeros only at boundaries,")+Zeros[l]+(Str("is inside interval (","")+lb+","+ub+")"));
            if(abs(lb-Zeros[l])<1.e-12)IZeros.push_back(Idx->basis()->integrable()->lowerMargin());
            if(abs(ub-Zeros[l])<1.e-12)IZeros.push_back(Idx->basis()->integrable()->upperMargin());
        }
    }
}

bool ConstrainedView::isBlockDiagonal() const{
    return _op->isBlockDiagonal();
}

bool ConstrainedZero::isDummy() const {
    if(_view) return dynamic_cast<OperatorDUM*>(oFloor) != nullptr;
    else return dynamic_cast<const OperatorDUM*>(dynamic_cast<Floor*>(oFloor)->_parentFloor) != nullptr;
}


const OperatorFloor* ConstrainedZero::floor() const {
    if(oFloor == nullptr or dynamic_cast<Floor*>(oFloor) == nullptr) return OperatorTree::floor();
    return dynamic_cast<const Floor*>(oFloor)->_parentFloor;
}

OperatorFloor*& ConstrainedZero::floor() {
    if(oFloor == nullptr or dynamic_cast<Floor*>(oFloor) == nullptr) return OperatorTree::floor();
    return dynamic_cast<Floor*>(oFloor)->_parentFloor;
}

ConstrainedZero::Floor::Floor(OperatorFloor *Parent, std::vector<unsigned int> LeftZeros, std::vector<unsigned int> &RightZeros)
    :OperatorFloor(Parent->rows(),Parent->cols(),"constrained["+Parent->kind()+"]"),
      _parentFloor(Parent),_leftZeros(LeftZeros),_rightZeros(RightZeros){}

ConstrainedZero::Floor::Floor(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf)
    : OperatorFloor("")
{
    std::size_t code = (std::size_t) Buf[0].real();
    auto info = Info;
    auto buf = Buf;
    info[0] = code;
    std::size_t numRightZeros = (std::size_t) Buf[1].real(),
                numLeftZeros  = (std::size_t) Buf[1].imag();
    _rightZeros.resize(numRightZeros);
    _leftZeros.resize(numLeftZeros);
    for(std::size_t n = 0; n < std::min(numLeftZeros, numRightZeros); ++n) {
        _rightZeros[n] = (std::size_t) Buf[2 + n].real();
        _leftZeros[n] = (std::size_t) Buf[2 + n].real();
    }
    if(numRightZeros > numLeftZeros) {
        for(std::size_t n = numLeftZeros; n < numRightZeros; ++n)
            _rightZeros[n] = (std::size_t) Buf[2 + n].real();
    }
    else if(numLeftZeros > numRightZeros) {
        for(std::size_t n = numRightZeros; n < numLeftZeros; ++n)
            _leftZeros[n] = (std::size_t) Buf[2 + n].imag();
    }

    auto prefix_size = 2 + std::max(numLeftZeros, numRightZeros);
    info[4] -= prefix_size;
    info[3] -= prefix_size;
    buf.erase(buf.begin(), buf.begin() + prefix_size);
    _parentFloor = unpackFactory(info, buf); // who is going to delete this memory?

    _rows = _parentFloor->rows();
    _cols = _parentFloor->cols();
    _kind = "constrained[" + _parentFloor->kind() + ']';
}

void ConstrainedZero::Floor::axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
                                  const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const {
    vector<complex<double> > xCopy;
    if(_rightZeros.size()>0){
        xCopy.assign(X,X+SizX);
        for(int k=0;k<_rightZeros.size();k++)xCopy[_rightZeros[k]]=0.;
        X=xCopy.data();
    }
    if(_leftZeros.size()==0)
        _parentFloor->apply(Alfa,X,SizX,Beta,Y,SizY);
    else {
        vector<complex<double> > rhs(SizY);
        _parentFloor->apply(Alfa,X,SizX,0.,rhs.data(),SizY);
        for(int k=0;k<_leftZeros.size();k++)rhs[_leftZeros[k]]=0.;
        for(int k=0;k<SizY;k++)Y[k]=Beta*Y[k]+rhs[k];
    }
}

void ConstrainedZero::Floor::pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const {
    Buf.emplace_back(packCode(_parentFloor->kind()), .0);
    Buf.emplace_back(_rightZeros.size(), _leftZeros.size());
    for(unsigned int k=0;k<std::min(_leftZeros.size(),_rightZeros.size());k++)Buf.push_back(complex<double>(_rightZeros[k],_leftZeros[k]));
    if(_leftZeros.size() > _rightZeros.size()) {
        for(unsigned int k = _rightZeros.size(); k < _leftZeros.size(); ++k)
            Buf.push_back(complex<double>(0.,_leftZeros[k]));
    }
    else if(_leftZeros.size() < _rightZeros.size()) {
        for(unsigned int k = _leftZeros.size(); k < _rightZeros.size(); ++k)
            Buf.push_back(complex<double>(_rightZeros[k],0.));
    }
    _parentFloor->pack(Info, Buf);
    Info[0] = packCode("constrained");
}

ConstrainedZeroDerivative::ConstrainedZeroDerivative(OperatorTree *Op, std::vector<std::string> Axes,
                                                     const std::vector<std::vector<double> > &Zeros)
{
    _op=Op;
    iIndex=Op->iIndex;
    jIndex=Op->jIndex;
    if(Op->isLeaf()){
        vector<vector<complex<double> > > uLeft,uRight;
        for(int k=0;k<Axes.size();k++){
            getTrans(Axes[k],Zeros[k],iIndex,uLeft);
            getTrans(Axes[k],Zeros[k],jIndex,uRight);
        }
        if((uLeft.size()>0 or uRight.size()>0))
            oFloor=new Floor(Op->floor(),uLeft,uRight);
        else
        {
            _view=true;
            oFloor=const_cast<OperatorFloor*>(Op->floor());
        }
    }
    for(int k=0;k<Op->childSize();k++)
        childAdd(new ConstrainedZeroDerivative(Op->child(k),Axes,Zeros));
}

void ConstrainedZeroDerivative::getTrans(std::string Axis, const std::vector<double> & Zeros, const Index* Idx,
                                         std::vector<std::vector<std::complex<double> > > &UTrans)
{
    if(Idx->heightAboveBottom()>0)ABORT("for now only 1-d constraints");

    if(Idx->axisName()==Axis and Idx->basis()->integrable()!=0){
        double lb=Idx->basis()->integrable()->lowBound();
        double ub=Idx->basis()->integrable()->upBound();
        for(int l=0;l<Zeros.size();l++){
            if(lb+1.e-12<Zeros[l] and Zeros[l]<ub-1.e-12)
                ABORT(Str("Zeros only at boundaries,")+Zeros[l]+(Str("is inside interval (","")+lb+","+ub+")"));
            if(abs(lb-Zeros[l])<1.e-12 or abs(ub-Zeros[l])<1.e-12){
                UTrans=Idx->basis()->integrable()->transZeroValDer(Zeros[l]);
                UTrans.resize(2);
            }
        }
    }
}

ConstrainedZeroDerivative::Floor::Floor(OperatorFloor *Parent,
                                        std::vector<std::vector<std::complex<double> > > &ULeft,
                                        std::vector<std::vector<std::complex<double> > > &URight)
    :OperatorFloor(Parent->rows(),Parent->cols(),"constrained["+Parent->kind()+"]"),
      _parentFloor(Parent),_uLeft(ULeft),_uRight(URight){}


void ConstrainedZeroDerivative::Floor::axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
                                            const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const {
    vector<complex<double> > xCopy;
    if(_uRight.size()>0){
        xCopy.assign(X,X+SizX);
        project(xCopy.data(),xCopy.size(),_uRight);
        X=xCopy.data();
    }
    if(_uLeft.size()==0)
        _parentFloor->apply(Alfa,X,SizX,Beta,Y,SizY);
    else {
        vector<complex<double> > lhs(SizY);
        _parentFloor->apply(Alfa,X,SizX,0.,lhs.data(),SizY);
        project(lhs.data(),lhs.size(),_uLeft);
        for(int k=0;k<SizY;k++)Y[k]=Beta*Y[k]+lhs[k];
    }
}

void ConstrainedZeroDerivative::Floor::project(std::complex<double>*Vec,int Size,
                                               const std::vector<std::vector<std::complex<double> > > & U) const
{
    for(int l=0;l<U.size();l++){
        std::complex<double> sum=0.;
        for(int k=0;k<Size;k++)sum+=std::conj(U[l][k])*Vec[k];
        for(int k=0;k<Size;k++)Vec[k]-=U[l][k]*sum;
    }
}

