// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef MULTIPOLEORBITAL_H
#define MULTIPOLEORBITAL_H

#include <complex>
#include <vector>
#include <string>
#include <map>
#include "qtEigenDense.h"
#include "coefficients.h"


class BasisOrbitalNumerical;
class Index;

/// organize BasisOrbital for use in multipole density
class MultipoleOrbital
{
    /// _listCoeff[Name][iChannel][m] ... vector of l-components
    static std::map<std::string,std::vector<std::map<int,std::vector<Coefficients*> > > > _listCoef;
    static std::map<std::string,std::vector<std::map<int,std::vector<int> > > > _listLangle;
    static std::map<std::string,Eigen::MatrixXcd> _listRho;

    const std::vector<Coefficients*>*_lCoef;
    const std::vector<int> *_lAngle;
    const Eigen::MatrixXcd  *_rho;
public:
    static int channelIndex(const std::string&Chan,const Index* Idx);
    static int lQuantumNumber(const Index* Idx);
    static int mQuantumNumber(const Index* Idx);
    /// organize coefficients by channel and m-quantum number
    static void add(std::string Name, const BasisOrbitalNumerical &Orbs,const std::vector<std::complex<double> > Rho={});
    /// all l-components for given channel and m-quantum number
    MultipoleOrbital(std::string Name, const Index* Idx);

    inline const Eigen::MatrixXcd vals(int K) const{
        return Eigen::Map<Eigen::MatrixXcd>(_lCoef->data()[K]->data(),
                                            _lCoef->data()[K]->size(),1);
    }
    int lSize() const{return _lAngle->size();}
    int lambda(int K) const{return _lAngle->data()[K];}
    int mMin() const;
    int mMax() const;
    int lMin() const;
    int lMax() const;
    std::complex<double> density(std::string Name, const Index* IIndex, const Index* JIndex) const; ///< matrix element of density matrix
};

#endif // MULTIPOLEORBITAL_H
