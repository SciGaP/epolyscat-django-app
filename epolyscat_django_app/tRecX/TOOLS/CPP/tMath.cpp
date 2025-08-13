// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "tMath.h"
#include "stdio.h"
#include <iostream>

long double tMath::doubleFactorial(unsigned int N){
    long double f=1.;
    if(N%2==0)for(unsigned int k=1;k<N/2+1;  k++)f*=(long double)(2*k);
    else      for(unsigned int k=1;k<(N+1)/2;k++)f*=(long double)(2*k+1);
    return f;
}

