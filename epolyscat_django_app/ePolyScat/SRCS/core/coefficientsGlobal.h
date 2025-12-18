// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COEFFICIENTSGLOBAL_H
#define COEFFICIENTSGLOBAL_H

#include "coefficients.h"
//#include "coefficientsFloor.h"

/** \ingroup Coefficients */

/// @brief Coefficients in storage sorted by Index ownership (data is NOT ordered in tree-ordering)
class CoefficientsGlobal: public Coefficients{
    static std::map<std::string,CoefficientsGlobal*> _views;
    std::vector<int> _sizes;
    std::complex<double>* _storageData;
    CoefficientsGlobal(const CoefficientsGlobal&){DEVABORT("not implemented");}
    void operator=(const Coefficients* NotAllowed);
public:
    static CoefficientsGlobal* view(Coefficients* C); ///< global view on Coefficients
    CoefficientsGlobal(Coefficients* C);
    CoefficientsGlobal(const Index* Idx, std::complex<double> Val=0.);
    std::vector<int> & sizes() {return _sizes;}
    CoefficientsGlobal& operator=(const Coefficients &);
    CoefficientsGlobal& operator=(const CoefficientsGlobal &);
    std::complex<double> * storageData(){return _storageData;}
    std::complex<double> * storageData() const {return _storageData;}
    std::string strNode(int Precision=0) const;

};

#endif // COEFFICIENTSGLOBAL_H
