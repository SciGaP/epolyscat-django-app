// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISCHANNEL_H
#define BASISCHANNEL_H

#include "qtEigenDense.h"
#include <memory>

class BasisOrbital;

#include "basisAbstract.h"

///@brief Channel basis Psi{I} defined through a set of orbitals and excitation density matrices rho{IJ}(i,j)
///
/// Psi{I} = sum[i1,...,in] A{I}[i1,...,iN] Phi[i1]*Phi[i2]*...*Phi[iN]
/// <br> Phi[i]...i'th BasisOrbital
/// <br> rho{IJ}(i,j) = sum[i2,...,iN] A{I}[i,i2,...iN] conj(A{J}[j,i2,...,iN])
class BasisChannel : public BasisAbstract
{
    std::vector<std::vector<Eigen::MatrixXcd>> _rho1;
    Eigen::ArrayXcd _rho2; // not in use for now
    const BasisOrbital*_orb;
public:
    unsigned int size() const {return _rho1.size();}
    BasisChannel(std::string Def /** ChanelXXX[OrbName:RefIdx], XXX=Det,HF,Hole */,
                 int Size, int First /** start hier in OrbName orbitals */);
    BasisChannel(std::string Def, const BasisOrbital * Orbs);
    const Eigen::MatrixXcd & rho(int I, int J) const {return _rho1[I][J];}
    const BasisOrbital * orbs() const {return _orb;}
    const std::vector<std::vector<Eigen::MatrixXcd>> & rho1() const {return _rho1;}
    void print() const;
    bool isOrthonormal() const {return true;}
    std::string strDefinition() const{return "undefined";}
};

#endif // BASISCHANNEL_H
