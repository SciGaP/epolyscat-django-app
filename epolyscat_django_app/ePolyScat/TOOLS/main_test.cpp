// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "readInput.h"
//#include "algebra.h"

#include <iostream>
#include <fstream>
#include "vectorReal.h"
#include "orthopol.h"
#include "orthogonalDerived.h"
#include "tools.h"
#include "odeTest.h"
#include "vectorComplex.h"
#include "treeExample.h"
#include "algebra.h"
#include "str.h"
#include "odeRungeKutta.h"

using namespace std;
using namespace tools;
using std::ofstream;
using std::ios_base;

int main(int argc, const char* argv[]) {

//     OrthogonalPolynomial::Test(true);
//    Algebra::Test();
//    Integrate::test();
//    Units::test();
//    ReadInput::Test(argc,argv);
//    TreeExample::test();
//    Str::Test();
//    OdeRungeKutta<TestExp,complex<double> >::test();

    OrthogonalPolynomial::Test(true);
    OrthogonalDerived::verify();



}
