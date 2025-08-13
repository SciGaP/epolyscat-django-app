// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMO_H
#define BASISMO_H

#include "basisAbstract.h"
#include <memory>
#include "vectorValuedFunction.h"

class ReadInput;
class mo;
class QuantumChemicalInput;

///@brief Molecular orbitals in cartesian coordinates
///
/// at present just a wrapper for Vinay Majety's old code
class BasisMO : public BasisAbstract, public VectorValuedFunction
{
    std::string _source;
    std::shared_ptr<mo> _vinayMO;
    std::vector<int> _select;
public:
    static void read(ReadInput & Inp);
    static void add(const QuantumChemicalInput * System);
    BasisMO(ReadInput &Inp);
    BasisMO(const QuantumChemicalInput * System);
    unsigned int size() const;
    unsigned int length() const{return size();}
    std::vector<std::complex<double> >operator()(std::vector<double> X)const;
    std::string coordinates() const {return "X.Y.Z";}
    void print(std::string RefIdx="main") const;
    std::string str(int Level=0) const;
    std::string info(int Number=-1) const;
    const mo* vinayMO() const {return _vinayMO.get();}
    bool isOrthonormal() const {return true;}
    std::string strDefinition() const{return "undefined";}
};

#endif // BASISMO_H
