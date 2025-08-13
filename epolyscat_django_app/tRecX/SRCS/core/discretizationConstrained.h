// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONCONSTRAINED_H
#define DISCRETIZATIONCONSTRAINED_H

#include "discretizationDerived.h"


/** \ingroup Discretizations */

/// reads input from file to constrain a given discretization
class DiscretizationConstrained: public DiscretizationDerived{
    std::string constString;
public:
    DiscretizationConstrained(const Discretization *D, ReadInput &Inp);
    std::string constraints() const {return constString;}

    static bool inputs(ReadInput & Inp); ///< parse all inputs
};

#endif // DISCRETIZATIONCONSTRAINED_H
