// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef MULTIPOLEPOTENTIAL_H
#define MULTIPOLEPOTENTIAL_H

#include <string>
#include <complex>
#include <memory>
#include <vector>
#include <map>
#include "qtEigenDense.h"
#include "abort.h"
#include "str.h"

class BasisIntegrable;

///@brief radial factors for the multipole expansion for V(|vec(r)-vec(s)|)
///
/// Matrix of exact integrals
///
///     V[l][i,j] = int dr ds conj(IBas[i](r)) V[l](r,s) JBas[j](s) * 4 pi/(2*l+1)
///
/// For V(r)=1/r (Coulomb):  V[l](r,s) = min(r,s)^l / max(r,s)^(l+1)
class MultipolePotential
{
public:
    typedef std::complex<double> (*multipoleRadial)(int Lambda, std::complex<double> R, std::complex<double> S);
private:
    static std::map<std::string,multipoleRadial> _listPots;
    std::shared_ptr<std::vector<Eigen::MatrixXcd> >_vals;
    static std::map<std::string,std::shared_ptr<std::vector<Eigen::MatrixXcd>>> _listVals;
    std::shared_ptr<std::vector<Eigen::MatrixXcd>> integrate(int LMax, multipoleRadial Pot, const BasisIntegrable *IBas, const BasisIntegrable *JBas);
public:
    MultipolePotential(int LMax /** for la=0...LMax */,
                       std::string PotRadial /** name-string to define potential */,
                       const BasisIntegrable * IBas, const BasisIntegrable * JBas);

    /// V[Lambda](i,j)
    const Eigen::MatrixXcd & vals(int Lambda) const
    {
        if(Lambda>=_vals->size())DEVABORT(Sstr+"multipole lambda exceeds stored size:"+Lambda+">="+_vals->size());
        return _vals->at(Lambda);
    }
    static void Test();
};

#endif // MULTIPOLEPOTENTIAL_H
