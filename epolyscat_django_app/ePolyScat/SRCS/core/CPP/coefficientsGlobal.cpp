// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "coefficientsGlobal.h"

#include "index.h"
#include "parallel.h"

#include "mpiWrapper.h"
#include "str.h"

using namespace std;

CoefficientsGlobal::CoefficientsGlobal(const Index* Idx, complex<double> Val):Coefficients(Idx){
    storageAssign(Idx->sizeStored(),Val);
    _storageData=Coefficients::storageData();

    // collect floors by floor owner
    vector<vector<Coefficients*> > fProc(MPIwrapper::Size());
    for(Coefficients * l=const_cast<Coefficients*>(firstLeaf());l!=0;l=l->nextLeaf()){
        if(Parallel::none==Parallel::owner(l->idx()))DEVABORT("no owner assigned");
        if(l->idx()->hasFloor())fProc[Parallel::owner(l->idx())].push_back(l);
    }

    // assign storage pointers
    unsigned int loc=0;
    _sizes.assign(MPIwrapper::Size(),0);
    for(unsigned int n=0;n<fProc.size();n++){
        for(unsigned int k=0;k<fProc[n].size();k++){
            fProc[n][k]->setFloorData(storageData()+loc);
            _sizes[n]+=fProc[n][k]->size();
            loc+=fProc[n][k]->size();
        }
    }
    unsetOrderedData(); // only floors are Index-ordered
}

/// create a Global view on C (storage of C will be rearranged)
CoefficientsGlobal::CoefficientsGlobal(Coefficients* C)
    :Coefficients(*C, [&](Coefficients* node){return node->isLeaf();})
{
    // collect floors by floor owner
    vector<vector<Coefficients*> > fProc(MPIwrapper::Size());
    for(Coefficients * l=const_cast<Coefficients*>(C->firstLeaf());l!=0;l=l->nextLeaf()){
        fProc[Parallel::owner(l->idx())].push_back(l);
    }
    // move content to Global contiguous and redirect pointers
    vector<complex<double> > stor(size());
    _sizes.assign(MPIwrapper::Size(),0);
    unsigned int loc=0;
    for(unsigned int n=0;n<fProc.size();n++){
        for(unsigned int k=0;k<fProc[n].size();k++){
            for(unsigned int l=0;l<fProc[n][k]->size();l++)
                stor[loc+l]=fProc[n][k]->floorData()[l];
            fProc[n][k]->setFloorData(stor.data()+loc);
            _sizes[n]+=fProc[n][k]->size();
            loc+=fProc[n][k]->size();
        }
    }
    C->replaceStorage(stor); // replace all C-storage with stor
    _storageData=C->storageData();

    // storage of C has been rearranged - not Index-ordered
    unsetOrderedData();
    C->unsetOrderedData();
}

CoefficientsGlobal& CoefficientsGlobal::operator=(const Coefficients & Other){
    if (this==&Other)return *this;
    if(orderedData()!=0 and Other.orderedData()!=0)
        memcpy(orderedData(),Other.orderedData(),size()*sizeof(*orderedData()));
    else
        for(size_t k=0;k<childSize();k++)child(k)->operator=(*Other.child(k));
    return *this;
}

CoefficientsGlobal& CoefficientsGlobal::operator=(const CoefficientsGlobal & Other){
    if (this==&Other)return *this;
    if(orderedData()!=0 and Other.orderedData()!=0)
        memcpy(orderedData(),Other.orderedData(),size()*sizeof(*orderedData()));
    else
        for(size_t k=0;k<childSize();k++)child(k)->operator=(*Other.child(k));
    return *this;
}

map<std::string,CoefficientsGlobal*> CoefficientsGlobal::_views;
CoefficientsGlobal* CoefficientsGlobal::view(Coefficients* C)
{
    if(_views.count(C->hash())==1)return _views[C->hash()];
    _views[C->hash()]=new CoefficientsGlobal(C);
    return _views[C->hash()];
}

string CoefficientsGlobal::strNode(int Precision) const {
    return Str("[","")+Parallel::owner(idx())+"] "+orderedData()+size();
}
