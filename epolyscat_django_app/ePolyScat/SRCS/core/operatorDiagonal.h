// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORDIAGONAL_H
#define OPERATORDIAGONAL_H

#include <vector>
#include <complex>
#include "abort.h"
#include "operatorAbstract.h"
#include "tree.h"

typedef std::complex<double> (*timeDependentFunction)(std::complex<double>,double);

/** \ingroup Structures */
/// diagonal operator and functions of it
class OperatorDiagonal : public Tree<OperatorDiagonal>,public OperatorAbstract
{
    std::vector<std::complex<double> > dVal;
    std::vector<std::complex<double> > funcVal;
    timeDependentFunction func;
    static std::complex<double> id(std::complex<double> Arg, double Time){return Arg;}
public:
    // a few frequent functions
    static std::complex<double> doNotUpdate(std::complex<double> Arg, double Time){ABORT("this should never be called");return Arg;}
    static std::complex<double> keepFunction(std::complex<double> Arg, double Time){ABORT("this should never be called");return Arg;}
    static std::complex<double> identityFunction(std::complex<double> Arg, double Time){return Arg;}
    static std::complex<double> zeroFunction(std::complex<double> Arg, double Time){return 0;}
    static std::complex<double> oneFunction(std::complex<double> Arg, double Time){return 1;}
    static std::complex<double> expItFunction(std::complex<double> Arg, double Time){return std::exp(std::complex<double>(0.,Time)*Arg);}

    ~OperatorDiagonal();
    OperatorDiagonal(){}
    OperatorDiagonal(const std::string &Name,const Index *Idx,timeDependentFunction Func=identityFunction);
    void apply(std::complex<double> Alfa, const Coefficients &X, std::complex<double> Beta, Coefficients &Y) const;

    void axpy(std::complex<double> Alfa, const Coefficients &X, std::complex<double> Beta, Coefficients &Y, const double Time);

    void add(const std::vector<std::complex<double> >& dVal,
             const Index * Idx,
             std::complex<double> Factor=1.);
    void update(double Time, const Coefficients* CurrentVec=0){updateFunction(Time,keepFunction);}
    void updateFunction(double Time, timeDependentFunction Func);
    void setFunctionValue(unsigned int K, std::complex<double> ValK); ///< directly set diagonal function values
    using Tree::str;
    std::string strNode(int Precision) const;
    const std::vector<std::complex<double> > & diagonal(){return dVal;}

    unsigned int vals() const; ///< number of diagonal values
    std::complex<double> val(unsigned int K) const; ///< the k'th eigenvalue
//    const std::vector<std::complex<double>> & allvals() const{return dVal;}; ///< the k'th eigenvalue
    void setVal(unsigned int K, std::complex<double> Val); ///< set the k'th eigenvalue

    /// Guarantees child(i)->iIndex == iIndex->child(i) for all children up until iIndex->hasFloor()
    void setupAccordingToIndex();
    /// Requires child(i)->iIndex == iIndex->child(i)
    void storeInCoefficients(Coefficients& c);
    /// Requires child(i)->iIndex == iIndex->child(i)
    void setFromCoefficients(const Coefficients& c);
};

#endif // OPERATORDIAGONAL_H
