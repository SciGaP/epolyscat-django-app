// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef JACOBIAN_H
#define JACOBIAN_H
#include "abort.h"
#include "useMatrix.h"
#include "str.h"
#include <memory>

class Jacobian{
public:
    virtual ~Jacobian(){}
    virtual void operator()(const UseMatrix & QuadX, UseMatrix & QuadW) const=0;
    virtual std::string kind() const=0;
    static const Jacobian* factory(std::string Kind,std::complex<double> Eta);
};

class Jacobian1:public Jacobian{
public:
    Jacobian1(std::complex<double> Eta){}
    std::string kind() const{return "1";}
    void operator()(const UseMatrix & QuadX, UseMatrix & QuadW) const{}
};

class JacobianQQ:public Jacobian{
public:
    std::string kind() const{return "QQ";}
    JacobianQQ(std::complex<double> Eta){}
    void operator()(const UseMatrix & QuadX, UseMatrix & QuadW) const{QuadW.cwiseProductIn(QuadX);QuadW.cwiseProductIn(QuadX);}
};

class JacobianQ:public Jacobian{
public:
    std::string kind() const{return "Q";}
    JacobianQ(std::complex<double> Eta){}
    void operator()(const UseMatrix & QuadX, UseMatrix & QuadW) const{QuadW.cwiseProductIn(QuadX);}
};

#endif // JACOBIAN_H
