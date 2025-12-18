// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisProd.h"

#include "index.h"
#include "indexGrid.h"
#include "vectorComplex.h"
#include "str.h"
//#include "basisMat.h"
#include "coordinateTrans.h"
#include "operatorDefinition.h"
#include "discretizationGrid.h"
#include "basisGrid.h"
#include "basisIntegrable.h"

std::map<const Index*,std::map<const Index*,BasisProd*>> BasisProd::all;

using namespace std;

Eigen::MatrixXd BasisProd::jacobianToNdim(const std::vector<double>&CoorNdim) const {
    return Eigen::Map<Eigen::MatrixXd>(_jacToNdim(CoorNdim).data(),dim(),dim());
}
std::vector<std::complex<double> > BasisProd::operator()(std::vector<double> CoorNdim) const{
    return BasisNdim::IndexAt(_idx.get(),CoorNdim);
}

BasisProd::BasisProd(const Index *RootIdx, const BasisNdim *Ref)
{
    _idx.reset(new Index(*RootIdx));
    _ndimCoor=RootIdx->coordinates();
    _fromNdim= CoordinateTrans::toCartesian(_ndimCoor);
    _toNdim=   CoordinateTrans::fromCartesian(_ndimCoor);
    _jacToNdim=CoordinateTrans::jacCartesian(_ndimCoor);
    _intFactor=CoordinateTrans::integrationFactor(_ndimCoor);
    _nabFactor=CoordinateTrans::nablaFactor(_ndimCoor);

    _quadCoor=Ref->quadCoor();
    CoordinateMap fromMain=CoordinateTrans::toCartesian(_quadCoor);
    IntegrationFactor mainFactor=CoordinateTrans::integrationFactor(_quadCoor);

    _valDer.resize(Ref->quadGrid().size());
    for(size_t pt=0;pt<Ref->quadGrid().size();pt++){
        // Ref._quadGrid,Ref._quadWeig is wrt to main coors, transform to Ndim coors
        _quadGrid.push_back(_toNdim(fromMain(Ref->quadGrid()[pt])));
        _quadWeig.push_back(Ref->quadWeig()[pt]*pow(mainFactor(Ref->quadGrid()[pt])/_intFactor(_quadGrid[pt]),2));

        // evaluate all basis functions at given point
        valDerAll(RootIdx,_quadGrid[pt],_valDer[pt]);
    }
    _comSca.resize(RootIdx->height()); // for now, no complex scaling on product basis

    // transform back to mainCoor
    mainQuadValDer();
}

BasisProd::BasisProd(const Index *RootIdx, unsigned int MinQuad)
{
    _idx.reset(new Index(*RootIdx));
    _ndimCoor=RootIdx->hierarchy_no_NONE();
    _quadCoor=_ndimCoor;
    productGridWeig(RootIdx,_quadGrid,_quadWeig,MinQuad);

    std::vector<std::vector<double>>weig,grid;
    for(const Index* rdx=RootIdx;rdx and rdx->axisName()!="NONE";rdx=rdx->descend()){
        grid.push_back({});
        weig.push_back({});
        size_t n=std::max(rdx->basis()->integrable()->order(),MinQuad);
        rdx->basis()->integrable()->quadRule(n,grid.back(),weig.back());
    }
    // a product grid
    _gdx.reset(new IndexGrid(RootIdx,tools::splitString(_quadCoor,'.'),grid,weig));

    _fromNdim= CoordinateTrans::toCartesian(_ndimCoor);
    _toNdim=   CoordinateTrans::fromCartesian(_ndimCoor);
    _jacToNdim=CoordinateTrans::jacCartesian(_ndimCoor);
    _intFactor=CoordinateTrans::integrationFactor(_ndimCoor);
    _nabFactor=CoordinateTrans::nablaFactor(_ndimCoor);

    _valDer.resize(_quadGrid.size());
    for(size_t pt=0;pt<_quadGrid.size();pt++){
        // evaluate all basis functions at given point
        valDerAll(RootIdx,_quadGrid[pt],_valDer[pt]);
    }
    if(!_valDer[0].size())DEVABORT("no values");
    _comSca.resize(RootIdx->height()); // for now, no complex scaling on product basis
}

const BasisProd* BasisProd::factory(const Index *Idx, const Index *Ref){
    if(Ref->basis()->ndim()==0)ABORT("must have BasisNdim on Ref'erence Index");

    // note: this check should be improved
    Idx=Ref->basis()->ndim()->rootIndex(Idx);
    if(Idx==0)ABORT("Index does not match Ref");

    if(all.count(Idx)==0 or all[Idx][Ref]==0){
        all[Idx][Ref]=new BasisProd(Idx,Ref->basis()->ndim());
    }
    return all[Idx][Ref];
}

vector<complex<double> > BasisProd::tensorMultiply(const vector<complex<double> > &LFac, const vector<complex<double> > &RFac) const {
    vector<complex<double> > res;
    for(size_t k=0;k<RFac.size();k++)
        for(size_t l=0;l<LFac.size();l++)
            res.push_back(LFac[l]*RFac[k]);
    return res;
}

void BasisProd::valDerAll(const Index *Idx, vector<double> Point, vector<vector<complex<double> > >&ValDer) const{
    if(Point.size()==0){
        ValDer.push_back(vector<complex<double> >(1,1.));
    }
    else {

        if(Idx->continuity()!=Index::npos){
            // for fem levels, insert index coordinate
            for(size_t k=0;k<Idx->childSize();k++){
                std::vector<std::vector<std::complex<double> > > valder;
                valDerAll(Idx->child(k),Point,valder);
                for(size_t l=0;l<valder.size();l++)ValDer.push_back(valder[l]);
            }
        }
        else {
            const BasisIntegrable* bas=Idx->basis()->integrable();
            if(not bas)DEVABORT("must have integrable basis, is: "+Idx->strNode());
            if(bas->complexScaling()->xScaled(Point[0]).imag()!=0.)
                ABORT("cannot have multi-coordinate basis (Ndim) in complex scaled region");
            UseMatrix uPt=UseMatrix::Constant(1,1,Point[0]),val,der;
            bas->valDer(uPt,val,der,true);
            //"zeroOutside" is broken - do by hand
            if(bas->lowBound()>Point[0] or Point[0]>bas->upBound()){
                val*=0;
                der*=0;
            }
            bool expand=false;
            if((expand=not bool(Idx->childSize())))Idx->bottomExpand();
            for(size_t k=0;k<Idx->childSize();k++){
                std::vector<std::vector<std::complex<double> > > valder;
                valDerAll(Idx->child(k),vector<double>(Point.begin()+1,Point.end()),valder);
                for(size_t l=0;l<valder.size();l++){
                    complex<double> tmp=valder[l][0];
                    for(size_t i=0;i<valder[l].size();i++)valder[l][i]*=val(k).complex();
                    // sort deriviatives such that highest level has lowest index (e.g. phi at valder[l][1]
                    valder[l].insert(valder[l].begin()+1,tmp*der(k).complex());
                    ValDer.push_back(valder[l]);
                }
            }
            if(expand)Idx->bottomUnexpand();

        }
    }
    if(ValDer.size()==0)DEVABORT("no point");
}
