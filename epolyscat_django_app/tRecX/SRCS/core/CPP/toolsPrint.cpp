// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "toolsPrint.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

#include "coefficients.h"
#include "index.h"

namespace ToolsPrint
{

void CoefByFloors(std::string FileName, const Coefficients *C){
    Coefficients c(C->idx());
    c=*C;
    c.purgeNearZeros(1.e-12);
    std::ofstream fil((FileName).c_str(),(std::ios_base::openmode) std::ios::beg);
    for(const Coefficients * f=c.firstLeaf();f!=0;f=f->nextLeaf(&c)){
        fil<<(Sstr+f->idx()->index()+std::vector<std::complex<double>>(f->anyData(),f->anyData()+f->size()))<<std::endl;
    }
}

}
