// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PRODUCTFUNCTION_H
#define PRODUCTFUNCTION_H

#include <memory>

#include "vectorValuedFunction.h"

class Algebra;
class ReadInput;

///@ingroup Functions
///@brief Vector-valued multi-argument function given by products of single-argument Algebra's
///
/// Format (example for radial polar coordinates Phi.Eta.Rn, H-atom)
/// <br> 1.1.Q*exp(-Q),1.Q.pow[2](Q)*exp(-Q/2),exp(i*Q).sqrt(1-Q*Q).pow[2](Q)*exp(-Q/2)...states m,l,n=(0,0,1),(0,1,2),(1,1,2)
class ProductFunction : public VectorValuedFunction
{
    std::string _coors;
    std::vector<std::vector<std::shared_ptr<Algebra> > >_algs;
public:
    static void read(ReadInput& Inp); /// read and add to list
    ProductFunction(std::string Coors /** e.g. X.Y.Z, Phi.Eta.Rn */,
                    std::vector<std::string> Algs /** list of dot-separated Algebra strings, one for each coordinate */);
    std::string coordinates() const {return _coors;}
    std::vector<std::complex<double> > operator()(const std::vector<double> X) const;
    unsigned int length() const{return _algs.size();}
};

#endif // PRODUCTFUNCTION_H
