// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COEFFICIENTSPERMUTE_H
#define COEFFICIENTSPERMUTE_H

#include "coefficients.h"

class Index;

/** \ingroup Coefficients */
/// \brief Coefficient for permuted index
///
/// accepts/returns data from/to Coefficient with original index sorting
class CoefficientsPermute : public Coefficients
{
    Coefficients _viewPermAsOrig; // view onto permuted coefficient with tree structure as original
    const Index *_indexOrig;
    void copyView(bool FromOrig, Coefficients *COrig, Coefficients *View, std::complex<double>*BData) const;
public:
    /// permuted coefficient will have the original level I at new level Perm[I]
    CoefficientsPermute(const Index* IOrig, std::vector<unsigned int> Perm=std::vector<unsigned int>());
    /// fill *this with permuted data from Corig (also return reference to *this)
    CoefficientsPermute & fromOrig(const Coefficients & COrig);
    /// fill COrig with back-permuted data from *this (also return reference to COrig)
    Coefficients & toOrig(Coefficients & COrig) const;
};

#endif // COEFFICIENTSPERMUTE_H
