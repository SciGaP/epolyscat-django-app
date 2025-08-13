// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMONOMIAL_H
#define BASISMONOMIAL_H

#include "basisIntegrable.h"
#include <vector>

///@brief monomial basis for tests
class BasisMonomial:public BasisIntegrable{
    int _order;
    int _from;
public:
    BasisMonomial(int Order, double Low, double Up, int From=0):BasisIntegrable(Low,Up),_order(Order),_from(From){_name="Monomial";}
    unsigned int size() const{return _order-_from+1;}
    unsigned int order() const{return _order;}
    void quadRule(int N, std::vector<double> & QuadX, std::vector<double> & QuadW) const{
        OrthogonalLegendre leg;
        std::vector<double> xx,ww;
        leg.quadrature(N,xx,ww);
        QuadX.clear();
        QuadW.clear();
        for(int k=0;k<N;k++){
            QuadX.push_back(lowBound()+(0.5*(1.+xx[k]))*(upBound()-lowBound()));
            QuadW.push_back(ww[k]*0.5*(upBound()-lowBound()));
        }
    }
    void valDer(const std::vector<std::complex<double> > & X,
                std::vector<std::complex<double> > & Val,
                std::vector<std::complex<double> > & Der, bool ZeroOutside=false) const{

        for(auto x: X){
            Val.push_back(std::pow(x,_from));
            Der.push_back(_from==0?0:double(_from)*std::pow(x,_from-1));
        }
        for(size_t k=_from+1;k<size();k++){
            for(size_t i=0;i<X.size();i++){
                Val.push_back(Val[X.size()*(k-1)+i]*X[i]);
                Der.push_back(Der[X.size()*(k-1)+i]*X[i]+Val[size()*(k-1)+i]);
            }
        }
    }
};

#endif // BASISMONOMIAL_H
