// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "histogram.h"

#include "qtEigenDense.h"

#include "readInput.h"
#include "index.h"
#include "coefficients.h"
#include "operatorDefinition.h"
#include "basisMat1D.h"
#include "asciiFile.h"
#include "coefficientsViewDeep.h"
#include "coorSystem.h"
#include "basisGrid.h"

using namespace std;

Histogram::Histogram(const Index *FromIndex, vector<string> InfStr, int Depth){

    if(InfStr.size()==0){
        string def;
        ReadInput::main.read("Histogram","definition",def,"",
                             "specify blank-separated list as 'FromCoor,ToCoor,BinCoor,Low,Up,Size FromCoor1,...",1,"histogram");
        InfStr=tools::splitString(def,' ');
    }

    if(InfStr.size()<=Depth)return;

    // get the current map
    Info info(InfStr[Depth]);
    if(nSibling()==0)_map=std::shared_ptr<Map>(new Map(FromIndex,info));
    else             _map=parent()->child(0)->_map;

    for(int k=0;k<info.size();k++){
        childAdd(new Histogram(FromIndex,InfStr,Depth+1));
        childBack()->_center=info._low+(k+0.5)*info._width;
    }
}

Histogram::Info::Info(string From_To_Which_Low_Up_Size){
    vector<string> in=tools::splitString(From_To_Which_Low_Up_Size,',');
    if(in.size()!=6)ABORT("need 6 comma-separated values, got: "+From_To_Which_Low_Up_Size);

    _from=shared_ptr<CoorSystem>(CoorSystem::factory(in[0]));
    _to  =shared_ptr<CoorSystem>(CoorSystem::factory(in[1]));
    if(_from->refSystem()!=_to->refSystem())ABORT("coordinates do not share reference coordinates: "
                                                  +_from->name()+"->"+_from->refSystem()+" "
                                                  +_to->name()+"->"+_to->refSystem()
                                                  );
    string s=in[1].substr(0,in[1].find(in[2]));
    _binI=std::count(s.begin(),s.end(),'.');
    _low=std::stod(in[3]);
    _size=std::stoi(in[5]);
    _width=(std::stod(in[4])-_low)/_size;
}
std::string Histogram::strNode(int Precision) const {
    if(not isLeaf()) return Tree::strNode(Tree_ptrsOnly);
    return Str("")+index()+":"+_value+"["+_center+"]";
}

string Histogram::Info::str() const {
    Str s("","");
    s+_from->name()+"->"+_to->name()+" bin["+_binI+"] on ["+_low+","+(_low+_size*+_width)+" n="+_size;
    return s;
}

Histogram::Map::Map(const Index* FromIndex, Info Hist, std::vector<double> Coor){
    if(not FromIndex->isLeaf()
            and FromIndex->basis()->grid()==0
            and (FromIndex->basis()->integrable()==0 or
                 not BasisMat1D("1",FromIndex->basis(),FromIndex->basis()).useMat().isIdentity(1.e-12)))
        ABORT("must be grid or orthonormal basis, is:"+FromIndex->basis()->str());

    if(FromIndex->isRoot()){
        // find the subset of index levels
        vector<string> subC=tools::splitString(Hist._from->name(),'.');
        vector<string> hier=tools::splitString(FromIndex->hierarchy(),'.');
        for(int k=0;k<subC.size();k++)
            _subLevel.push_back(find(hier.begin(),hier.end(),subC[k])-hier.begin());
    }

    if(Coor.size()==Hist._from->dim()){

        std::vector<double> toCoor=Hist._to->fromRef(Hist._from->toRef(Coor));

        // compute weight
        double weight=Hist._to->jacRefdCoor(toCoor).determinant()
                /     Hist._from->jacRefdCoor(Coor).determinant();

        // assign histogram bin (>size() indicates outside)
        int bin=int((toCoor[Hist._binI]-Hist._low)/Hist._width);
        if(bin<0)bin=Hist.size();
        //NOTE: need to clarify whether inverse should be used
        _dat=std::shared_ptr<Data>(new Data(bin,abs(1./weight),Coor));
        return;
    }

    if(FromIndex->isLeaf())
        ABORT("coordinates "+Hist._from->name()+" not in hierachy "+FromIndex->root()->hierarchy());

    while(Hist._from->name().find(FromIndex->axisName())==string::npos)FromIndex=FromIndex->descend();
    if(FromIndex->basis()->grid()==0)ABORT("FromCoor axis not grid, basis: "+FromIndex->basis()->str());

    Coor.push_back(0.);
    for (int k=0;k<FromIndex->childSize();k++){
        Coor.back()=FromIndex->basis()->grid()->mesh()[k];
        childAdd(new Map(FromIndex->child(k),Hist,Coor));
    }
    Coor.pop_back();
}
string Histogram::Map::Data::str() const {
    return Str("","")+_bin+": weight="+_weight+", coor: "+_fromCoor;
}
string Histogram::Map::strNode(int Level) const {
    if(not isLeaf())return strNode(Tree_ptrsOnly);
    return _dat->str();
}

void Histogram::fill(Coefficients *GridVals){
    if(empty())return;

    // create deep view
    CoefficientsViewDeep view(GridVals->idx());
    for(Histogram * h=firstLeaf();h!=0;h=h->nextLeaf())h->_value=0.;

    fill(view.view(GridVals),std::vector<unsigned int>());
}

void Histogram::fill(const Coefficients * GridVals, std::vector<unsigned int> Idx){
    if(GridVals->isLeaf()){
        add(Idx,std::norm(GridVals->data()[0]));
    }
    else {
        Idx.push_back(0);
        for (int k=0;k<GridVals->childSize();k++){
            Idx.back()=k;
            fill(GridVals->child(k),Idx);
        }
    }
}

void Histogram::add(const std::vector<unsigned int> Idx, double Val){
    std::vector<unsigned int> sub;
    if(isLeaf()){
        _value+=Val;
    }
    else {
        // current sub-index
        for(unsigned int k=0;k<_map->_subLevel.size();k++)
            sub.push_back(Idx[_map->_subLevel[k]]);

        // select bin according to subset of Idx
        int b;
        if((b=bin(sub))<childSize())
            child(b)->add(Idx,Val*weight(sub));
    }
}

void Histogram::write(string FileName){
    // create AsciiFile
    AsciiFile afil(FileName);

    // write header
    afil.writeComments({"#  sum[Phi,Eta]: (specRn) = (0)"});

    std::vector<double> row(height()+1);
    write(afil,row);

}

void Histogram::write(AsciiFile & File, std::vector<double> & Row){
    if(height()==0 and depth()==0)return;
    if(depth()>0)Row[depth()-1]=_center;
    if(isLeaf()){
        Row.back()=_value;
        File.writeRow(Row);
    }
    for(int k=0;k<childSize();k++){
        child(k)->write(File,Row);
        //        if(depth()>1)File.writeBlankRow();
    }
}
