// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISSETDEF_H
#define BASISSETDEF_H

#include "coordinate.h"
#include "complexScaling.h"
#include <string>

class Index;
class BasisSetDef{
    friend class IndexNew;
    friend class Axis;
    std::string _modify; // instructions for path-dependent modifications
public:
    unsigned int order;
    unsigned int minOrder; // use functions starting from minOrder
    double shift,scale;
    std::string funcs;
    bool first,last;
    Coordinate coor;
    std::vector<int> margin;
    ComplexScaling comSca;
    bool exactIntegral; // perform integrations exactly
    bool deriv;
    std::vector<double> par;

    virtual ~BasisSetDef(){}
    BasisSetDef():order(0),minOrder(0),shift(0.),scale(1.),funcs("vector"),first(false),last(false),coor(),comSca(ComplexScaling()),exactIntegral(false),deriv(false){}
    BasisSetDef(unsigned int Order, double Shift, double Scale,std::string Funcs,bool ExactIntegral,
                const bool First, const bool Last, const Coordinate & Coor, std::vector<int> Margin={0,-1},
                ComplexScaling ComSca = ComplexScaling(), bool Deriv = false, std::vector<double> Par=std::vector<double>(0))
        : _modify(""),order(Order),minOrder(0),shift(Shift),scale(Scale),funcs(Funcs),first(First),last(Last),coor(Coor),margin(Margin),
          comSca(ComSca),exactIntegral(ExactIntegral),deriv(Deriv),par(Par){
        if(coor.isPeriodic()){
            first=false;
            last=false;
        }

        // modify parameters: interprete first element as reverted, if it has infinite range
        if(first and not last and BasisFunction::asympZero(funcs) and scale>0){
            shift+=scale;
            scale=-std::abs(scale);
        }
        // no shift or scale for useIndex
        if(funcs=="vector"){
            shift=0.;
            scale=1.;
        }
        if(margin[1]==-1)margin[1]=order-1;
    }
    /// resolve dependencies of basis definition on other functions
    BasisSetDef resolveDependence(const std::vector<unsigned int> &Pos, const std::vector<const Index *> &Path) const;

    std::string size() const;

    std::string name() const {return funcs;}
//    void resolveDependence(std::vector<const BasisSet*> &Bas, const std::vector<unsigned int> &N);

    double lowBound() const;
    double upBound() const;
    bool operator==(const BasisSetDef &o) const;
    std::string str() const;
    std::string modify() const{return _modify;};
};
#endif // BASISSETDEF_H
