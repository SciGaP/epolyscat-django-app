// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATOR_ABSTRACT_PRODUCT_H
#define OPERATOR_ABSTRACT_PRODUCT_H

#include <string>
#include <vector>
#include <complex>

#include "operatorAbstract.h"
#include "coefficients.h"

class OperatorTree;


class OperatorAbstractProduct: public OperatorAbstract{
private:
    std::vector<Coefficients> coeffs;

public:
    std::vector<const OperatorAbstract*> ops;
    
    OperatorAbstractProduct(std::string name, std::vector<const OperatorAbstract*> Ops);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

    long applyCount() const;
};





#endif //OPERATOR_ABSTRACT_PRODUCT_H
