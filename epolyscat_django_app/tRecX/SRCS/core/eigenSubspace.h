// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef EIGENSUBSPACE_H
#define EIGENSUBSPACE_H

#include <complex>
#include <vector>
#include <coefficientsMulti.h>


//class Operator;
class Discretization;
class DiscretizationSpectral;
class ReadInput;

class EigenSubspace
{
    double tolerance;

    const OperatorAbstract* h;
    DiscretizationSpectral * specH0;

    std::vector<std::vector<std::complex<double> > > e0;
    CoefficientsMulti c0;

    std::vector<std::complex<double> > eGuess;  ///< guess for energies including interaction
    std::vector<std::complex<double> > eResolv; ///< energies to be removed from resolvent

    std::complex<double> eCenter() const; ///< average of eResolv
    bool converged(const UseMatrix &Smat);
    CoefficientsMulti & applyResolvent(CoefficientsMulti & C) const;
    /// Sort - sorting of Values by their smallest distance to any of the Centers
    void sortByDistance(const std::vector<std::complex<double> > &Centers, const UseMatrix & OTarget, const std::vector<std::complex<double> > &Values, std::vector<int> & Sort);
    /// prepare precondintioner for given ETarget energies, with E0 resolvent energies
    void setTarget(std::vector<std::complex<double> > E0, std::vector<std::complex<double> > &ETarget,  std::vector<Coefficients *> & Evec);
    /// run to convergence
    std::vector<std::complex<double> > iterate();
public:
    ~EigenSubspace(){}
    EigenSubspace(ReadInput & Inp);
    void setup(const Discretization *D, const OperatorAbstract *H, const OperatorAbstract *H0);
    void eigen(const std::vector<std::complex<double> > & E0, std::vector<std::complex<double> > &ETarget, std::vector<Coefficients*> & Evec);
    void write(const std::string & File);
    void print() const;
};




#endif // EIGENSUBSPACE_H
