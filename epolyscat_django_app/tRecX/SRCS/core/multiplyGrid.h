// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef MULTIPLYGRID_H
#define MULTIPLYGRID_H

#include <string>
#include <complex>
#include <vector>
#include <cfloat>
#include "tree.h"
//#include "operator.h"
#include "operatorAbstract.h"

class Discretization;
class DiscretizationDerived;
class CoefficientsFunction;
class Index;

/// \ingroup Structures
///@brief Apply multiplication operator by transforming to quadrature grid and back
class MultiplyGrid:public OperatorAbstract
{
public:
    ~MultiplyGrid(){}
    MultiplyGrid();
    MultiplyGrid(const Discretization *D /** original, non-grid discretization */,
                 std::string Func /** string specifying the function (see "function factory" in code) */,
                 std::string Grid /** axis to transform: format as axisName1.axisName2.axisName3 */);

    void apply(std::complex<double> Alfa, const Coefficients & X, std::complex<double> Beta, Coefficients & Y) const;
    double setNorm(){return 1.;}
private:
    DiscretizationDerived* dGrid;
    CoefficientsFunction* function;

};

#endif // MULTIPLYGRID_H
