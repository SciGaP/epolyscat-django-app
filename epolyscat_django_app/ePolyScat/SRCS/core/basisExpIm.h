// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISEXPIM_H
#define BASISEXPIM_H
class BasisSub;

#include "basisIntegrable.h"

class BasisTrigon : public BasisIntegrable
{
protected:
    std::vector<int> _mValues;
    BasisTrigon(const BasisSetDef & Def);
    BasisTrigon(std::string StrDefinition);
    BasisTrigon(std::vector<int> MValues);
public:

    /// exp(i m ...) or m>0 sin(m..), else cos(m...)
    double physical(int N) const {return double(_mValues[N]);}
    int mValueAbs(int N) const {return std::abs(_mValues[N]);}

    unsigned int size() const {return _mValues.size();}
    unsigned int order() const;
    std::string str(int Level) const;
    std::string strDefinition() const;

    void valDer(const std::vector<std::complex<double> > & X,
                std::vector<std::complex<double> > & Val,
                std::vector<std::complex<double> > & Der, bool ZeroOutside=false) const=0;
    void quadRule(int N, std::vector<double> & QuadX, std::vector<double> & QuadW) const;
    bool operator==(const BasisAbstract& Other) const;
};

class BasisExpIm: public BasisTrigon {
public:
    BasisExpIm(const BasisSetDef & Def):BasisTrigon(Def){}
    BasisExpIm(std::string StrDefinition):BasisTrigon(StrDefinition){_name="ExpIm";}

    int mValue(int N) const {return _mValues[N];}
    void valDer(const std::vector<std::complex<double> > & X,
                std::vector<std::complex<double> > & Val,
                std::vector<std::complex<double> > & Der, bool ZeroOutside=false) const;
};

class BasisCosSin: public BasisTrigon {
public:
    BasisCosSin(const BasisSetDef & Def):BasisTrigon(Def){}
    BasisCosSin(std::string StrDefinition):BasisTrigon(StrDefinition){_name="CosSin";}
    void valDer(const std::vector<std::complex<double> > & X,
                std::vector<std::complex<double> > & Val,
                std::vector<std::complex<double> > & Der, bool ZeroOutside=false) const;
};

#endif // BASISEXPIM_H
