// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONFACTOR_H
#define DISCRETIZATIONFACTOR_H

#include <vector>
#include "discretizationDerived.h"
#include "basicDisc.h"

class ReadInput;

/** \ingroup Discretizations */
/// Factor discretization extracted from basic disc
class DiscretizationFactor : public BasicDisc
{

    enum FactorizationType{
      Ion_Residual,
      Residual_Ion
    };
    FactorizationType factorization;

    bool hasChannelLevel;    // With or without channel level

    /// Helpers for contract, if channel disc, i.e complement=true
    Coefficients *tempC1,*helper_ion;
    Coefficients *tempParentC, *tempViewC;
    Index* resIdx;    
    void initialize_helpers(const Discretization *Disc);

    std::complex<double> computeInnerProduct(const Coefficients & C, const Coefficients & Ion);          ///< Helper in the case of axis dependencies

    void contractEachChannel(const Coefficients & C, const Coefficients & Ion, Coefficients & Cofactor); ///< contract C with FactorC

public:    
    DiscretizationFactor(const Discretization * Disc, std::string Axes, bool Complement=false, int noChannels=0);
    ~DiscretizationFactor();
    void contract(const Coefficients & C, const Coefficients & FactorC, Coefficients & Cofactor);        ///< contract C with FactorC
};

#endif // DISCRETIZATIONFACTOR_H
