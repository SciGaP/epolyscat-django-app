// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORLOWRANK_H
#define OPERATORLOWRANK_H

#include "operatorSingle.h"
class CoefficientsFloor;

/** \ingroup Structures */
/// (incomplete)
class OperatorLowRank : public OperatorSingle { //i think this object does not work as intended, and it is used nowhere used
public:
    OperatorLowRank(const std::string & name, const std::string & definition, Discretization * iDisc, Discretization * jDisc,
                    std::vector<int> & iBlock, std::vector<int> & jBlock);
    ~OperatorLowRank();
    virtual void axpy(CoefficientsFloor & X,  CoefficientsFloor & Y, std::complex< double > & factor) const;
    void inverse();
    bool isZero(double Eps=0.) const;
};


#endif // OPERATORLOWRANK_H
