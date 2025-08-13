// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISEIGEN_H
#define BASISEIGEN_H

#include "basisOrbital.h"
#include "basisSetDef.h"
class VectorValuedFunction;
class Index;

/** \ingroup Basissets */

///@brief BasisOrbital composed of eigenvectors of a given operator
///
/// name has the format Eigenbasis[operator:reference], where
///
/// operator... convertible by OperatorDefinitionNew<br>
/// reference.. name of a refererence Index subset that has been set up before
///
/// Eigenvalues are sorted by increasing real part and the NLow'th eigenvector is the first orbital
class BasisEigen : public BasisOrbital
{
    std::string _operDef;
    std::string _refName;
    std::vector<int> _select;
    std::vector<std::complex<double> > _eigenValues;
    BasisEigen(const std::string Def);
public:
    static BasisEigen *factory(std::string Def);
    BasisEigen(std::string OperRef, int Size, int NLow=0);///<OBSOLESCENT
    BasisEigen(const BasisSetDef & Def);
    void generateOrbitals(const Index* Idx=0);
    std::string strDefinition() const;
    const std::vector<std::complex<double>>  & eigenvalues() const{return _eigenValues;}
};

#endif // BASISEIGEN_H
