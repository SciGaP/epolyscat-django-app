// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorAbstractProduct.h"

#include "index.h"
#include "operatorTree.h"


OperatorAbstractProduct::OperatorAbstractProduct(std::string name, std::vector<const OperatorAbstract*> Ops):
    OperatorAbstract("Product("+name+")", Ops[0]->iIndex, Ops[Ops.size()-1]->jIndex), ops(Ops){

    for(int i=1; i<ops.size(); i++){
        if(ops[i-1]->jIndex!=ops[i]->iIndex) ABORT("Index mismatch");
        coeffs.push_back(Coefficients(ops[i]->iIndex));
        coeffs.back().treeOrderStorage();
    }
}
    

void OperatorAbstractProduct::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    bool zero=false;

    ops[ops.size()-1]->apply(1., Vec, 0., const_cast<Coefficients&>(coeffs[coeffs.size()-1]));
    for(int i=ops.size()-2; i>0; i--){

        // Prevent application of further operators if the result is zero
        if(coeffs[i].isZero()){
            zero=true;
            break;
        }

        ops[i]->apply(1., coeffs[i], 0., const_cast<Coefficients&>(coeffs[i-1]));
    }

    if(zero){
        ops[0]->apply(0, coeffs[0], B, Y);
    }else{
        ops[0]->apply(A, coeffs[0], B, Y);
    }

}

long OperatorAbstractProduct::applyCount() const{
    long result=0;
    for(int i=0; i<ops.size(); i++) result+=ops[i]->applyCount();
    return result;
}
