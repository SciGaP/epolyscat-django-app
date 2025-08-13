// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ORTHOGONALDERIVED_H
#define ORTHOGONALDERIVED_H

#include <memory>
#include "orthopol.h"

class Algebra;
class OrthogonalDerived: public OrthogonalPolynomial{
    std::vector<std::vector<double> >recCoe;
    std::vector<double> _normsq;
    std::shared_ptr<OrthogonalPolynomial>_base;
    std::shared_ptr<Algebra> _aWeig,_aDer;

public:
    OrthogonalDerived(unsigned int MaxDegree, const OrthogonalPolynomial *Base, std::string Weig, std::string Der);
    inline double normsq(int I) const{return _normsq[I];}
    inline long double lowerBoundary() const {return _base->lowerBoundary();}
    inline long double upperBoundary() const {return _base->upperBoundary();}
    double weight(double X) const;
    double derWeight(double X) const;

    static void verify(); // test functionality
    typedef double (*positiveFunction)(double);
private:
    void construct(std::vector<std::vector<double> > Ovr);
    inline long double a(int i) const {return recCoe[0][i-1];}
    inline long double b(int i) const {return recCoe[1][i-1];}
    inline long double c(int i) const {return recCoe[2][i-1];}

    std::vector<std::vector<double> > overlap(unsigned int Degree);

};


#endif // ORTHOGONALDERIVED_H
