// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "indexCoulX.h"
#include "basisBesselCoulomb.h"
#include "operatorDefinition.h"
#include "inverseCoulX.h"
#include "bandedOverlap.h"
#include "basisAbstract.h"

IndexCoulX::IndexCoulX(const Index* Idx, double RC, double RMax, const std::vector<double> KGrid, bool BandOvr, bool PureBessel){
    setBasis(Idx->basis());
    setAxisName(Idx->axisName());
    if(Idx->isRoot()==true){
        fromIndex.push_back(Idx);
    }
    // needed for "surfaceDiscCoulX"
    if(Idx->axisName().find("surfRn")==0){
        setAxisName("Rn");
    }

    if(Idx->axisName().find("ValDer")==0){
        unsigned int lAngular=0;
        for(const Index *idxPar = Idx;idxPar!=Idx->root();idxPar=idxPar->parent()){
            if(idxPar->parent()->axisName().find("Eta")==0){
                lAngular = idxPar->parent()->basis()->physical(idxPar->nSibling());
                break;
            }
        }
        setBasis(new BasisBesselCoulomb(RC,RMax,lAngular,KGrid));

        setAxisName("Rn");
        setFloor(0);
        for(unsigned int l=0;l<KGrid.size();l++){
            childAdd(new Index(*Idx->child(0)));
        }
    }

    else{
        for(unsigned int k=0;k<Idx->childSize();k++){
            childAdd(new IndexCoulX(Idx->child(k),RC,RMax,KGrid,BandOvr,PureBessel));
        }
    }

    sizeCompute();
    setOverlap(new OperatorTree("ovr",OperatorDefinition("<<1>>",hierarchy()),this,this));

    unsigned int subD(0), superD(0);

    if(BandOvr){ // use banded overlap
        subD = 5; superD = 5; // Which values should they have?
        setOverlap(new BandedOverlap(this,subD,superD));
    }

    setInverseOverlap(new InverseCoulX(this,subD,superD,BandOvr));
}
