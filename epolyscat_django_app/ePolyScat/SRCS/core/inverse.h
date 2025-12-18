// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INVERSE_H
#define INVERSE_H

#include "operatorTree.h"
#include "coefficientsLocal.h"

class IIndex;
class JIndex;
/// Inverse composed of direct inverse S0^-1 and correction: S^-1(v) = correction( S0^-1(v) )
class Inverse: public OperatorTree
{
protected:
//    const OperatorAbstract* s0inv;
public:
    static const Inverse *factory(const Index *Idx);
//    Inverse(const OperatorAbstract* S0inv):OperatorTree(S0inv->name,S0inv->iIndex,S0inv->jIndex),s0inv(S0inv){}
//    Inverse(std::string Name, const Index* IIndex, const Index* JIndex):OperatorTree(Name,IIndex,JIndex),s0inv(0){}
    Inverse(const OperatorAbstract* S0inv):OperatorTree(S0inv->name,S0inv->iIndex,S0inv->jIndex){}
    Inverse(std::string Name, const Index* IIndex, const Index* JIndex):OperatorTree(Name,IIndex,JIndex){}

    virtual void applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const=0;
    virtual void apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const=0;
    virtual void project(Coefficients & Vec) const {} ///< project onto vector where inverse is defined

    /// for parallel application: use CoefficientsLocal
    virtual void parallelSetup() const {ABORT("no parallel setup implmented");}
    virtual void applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const=0;
    virtual void apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const=0;

    virtual void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    virtual void apply(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const;

    void verify(const OperatorAbstract * Ovr) const; ///< print matrices S Sinv and Sinv S
    void assignToIndex(Index* Idx) const;
//    void update(double Time, const Coefficients* C=nullptr) override {}
};

#endif // INVERSE_H
