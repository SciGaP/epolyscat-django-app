// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COEFFICIENTSLOCAL_H
#define COEFFICIENTSLOCAL_H
#include "coefficients.h"

#include <map>
#include <memory>

class Index;
class CoefficientsGlobal;
class ParallelContinuity;

/** \ingroup Coefficients */

/// @brief only data that are local on a given thread are stored
///
/// contiguous  and ordered storage matches the corresponding section in CoefficientsGlobal
class CoefficientsLocal: public Coefficients{
    static std::map<std::string,CoefficientsLocal*> _views;
    static std::complex<double> _dummyStorage;
    long unsigned int _size;
    static int setSize(const Coefficients* C);
public:
    std::shared_ptr<ParallelContinuity> newCont;
public:
    static bool local(const Coefficients*C); ///< true if C->idx() is on local process or not owned by any process
    static CoefficientsLocal* view(Coefficients* C); ///< local view on Coefficients
//    static CoefficientsLocal* view(CoefficientsGlobal* C); ///< local view on Global (keep data in place, i.e. Global remains intact)

    CoefficientsLocal(const CoefficientsLocal & Other);
    CoefficientsLocal():_size(0){}
    CoefficientsLocal(const Index *I, std::complex<double> Val=0.); ///< local Coefficient
    std::complex<double> * storageData();
    long unsigned int size() const override {return _size;}
    double norm() const override;
    std::string strNode(int Precision) const override;

    double localNorm() const {return Coefficients::norm();}
    CoefficientsLocal & operator=(const CoefficientsLocal & Other);
    void makeContinuous(double Scal=1.) override;
    CoefficientsLocal * cwiseProduct(const CoefficientsLocal & Rhs);


};

#endif // COEFFICIENTSLOCAL_H
