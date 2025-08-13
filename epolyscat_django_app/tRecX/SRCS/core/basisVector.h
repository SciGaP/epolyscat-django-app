// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISVECTOR_H
#define BASISVECTOR_H

#include "basisAbstract.h"
#include <map>
#include <memory>
#include "stringTools.h"

/** \ingroup Basissets */

///@brief Vector of complex numbers, replaces BasisSet(function useIndex)
class BasisVector : public BasisAbstract
{
    static std::map<std::string,std::unique_ptr<const BasisVector> >_list;
    unsigned int _size;
public:
    BasisVector(int Size):BasisAbstract("Vector"),_size(Size){}
    unsigned int size() const override{return _size;}
    static const BasisVector* factory(std::string Def){
        auto pBas=_list.find(Def);
        if(pBas==_list.end()){
            if(Def.find("Vector:")!=0)ABORT("need format Vector:size, got: "+Def);
            _list[Def]=std::unique_ptr<const BasisVector>(new BasisVector(tools::string_to_int(Def.substr(Def.find(":")+1))));
            pBas=_list.find(Def);
        }
        return pBas->second.get();
    }
    std::string strDefinition() const override {return "Vector:"+tools::str(size());}
    bool operator==(const BasisAbstract &other) const override {
        const BasisVector* o=dynamic_cast<const BasisVector*>(&other);
        return o?o->_size==_size:false;
    }
    bool isIndex() const override {return true;}
    bool isOrthonormal() const override {return true;}
    std::string str(int Level=0) const override {return strDefinition();}
};

#endif // BASISVECTOR_H
