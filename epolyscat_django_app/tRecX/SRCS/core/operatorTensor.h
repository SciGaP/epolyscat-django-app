// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORTENSOR_H
#define OPERATORTENSOR_H

#include "operatorSingle.h"

class OperatorTensor : public OperatorSingle {
public:
    /// TOO MANY CONSTRUCTORS
    OperatorTensor(){}
    OperatorTensor(const std::string & Def, const Index *IFloor, const Index *JFloor, std::vector<UseMatrix> & mats);
    OperatorTensor(const std::string & Def, const Discretization * iDisc, const Discretization * jDisc, const Index * IFloor, const Index * JFloor,
                   std::complex<double> Multiplier=1.);

    virtual ~OperatorTensor();

    /// Y = Alfa * Operator*X + Beta * Y;
    //    void axpy(const std::complex<double> Alfa,CoefficientsFloor & X, const std::complex<double> Beta, CoefficientsFloor & Y, bool transpose=false) const;
    void axpy(const std::complex<double> Alfa,CoefficientsFloor & X, const std::complex<double> Beta, CoefficientsFloor & Y, bool transpose=false) const{
        DEVABORT("not implemented");
    }
    void axpy(CoefficientsFloor & X,  CoefficientsFloor & Y, bool transpose=false) const{axpy(1.,X,1.,Y,transpose);} ///< Y = Operator*X + Y;
    void apply(std::complex<double> * InOut) const; ///< Inout = Operator*InOut;
    void apply(std::vector<std::complex<double> > & InOut) const; ///< Inout = Operator*InOut;
    void apply(CoefficientsFloor * X, CoefficientsFloor * Y, std::complex<double> Alfa=1.,std::complex<double> Beta=1.) const; ///< Y = alfa Op X + beta Y
    void inverse(); ///< replace all tensor factors by their inverses
    bool isZero(double Eps=0.) const; ///< true if OperatorTensor is zero
    std::string str() const; ///< return string describing the object
    std::string strDataStructure() const; ///< brief data structure
    void matrix(UseMatrix & Mat) const;
    void setStructureFlags(const std::string & Definition);

    //    virtual void packInfo(std::string &Info);
    //    virtual void unpackInfo(std::string &Info);

protected:
    enum matrixStructure {Id,Bd,Fu,IdFu,FuId,FuFu};
    matrixStructure mStruc;
    void initialize();
    bool absorb(OperatorSingle *&Other);
};


#endif // OPERATORTENSOR_H

