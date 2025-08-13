// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "subIter.h"

SubIter::SubIter()
{
}


void SubIter::eigen(unsigned int Nvec, std::vector<Coefficients *> Rvec, std::complex<double> Eval, double Eps){
    // check Rvec, possibly fill with random values and orthogonalize

    // extend to working size Nwork

    // iterate:
    // generate new set of vectors from old by applying preconditioned operator

    // solve larger eigenproblem

    // reduce to Nwork

    // convergence criterion: check residues
}
