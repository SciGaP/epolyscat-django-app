// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISASSOCLEG_H
#define BASISASSOCLEG_H

#include "basisIntegrable.h"

class BasisAssocLeg: public BasisIntegrable
{
    int _mAbs;
    int _size;
public:
    BasisAssocLeg(std::string Def);
    std::string name() const{return "assocLegendre["+tools::str(_mAbs)+"]";}
    unsigned int size() const {return _size;}
    unsigned int order() const {return _size+_mAbs+(_size+_mAbs)%2+1;}
    double physical(int Index) const {return _mAbs+Index;}
    double lValue(int Index) const {return _mAbs+Index;}
    double lowBound() const{return -1.;}
    double upBound() const{return 1.;}

    std::string str(int Level=0) const {return name()+" ["+tools::str(size())+"]";}
    static std::string strDefinition(const BasisSetDef& Def);
    std::string strDefinition() const {return std::string("AssociatedLegendre: ")+tools::str(_mAbs)+","+tools::str(size());}

    void valDer(const std::vector<std::complex<double> > & X,
                std::vector<std::complex<double> > & Val,
                std::vector<std::complex<double> > & Der, bool ZeroOutside=false) const;
    void quadRule(int N, std::vector<double> & QuadX, std::vector<double> & QuadW) const;
    bool operator==(const BasisAbstract& Other) const {return strDefinition()==Other.strDefinition();}
};

#endif // BASISASSOCLEG_H
