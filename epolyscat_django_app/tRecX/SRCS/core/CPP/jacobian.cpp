// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "jacobian.h"

const Jacobian* Jacobian::factory(std::string Kind,std::complex<double> Eta)
{
    if(Kind=="1" )return new Jacobian1(Eta);
    if(Kind=="Q" )return new JacobianQ(Eta);
    if(Kind=="QQ")return new JacobianQQ(Eta);
    DEVABORT("undefined Jacobian string: \""+Kind+"\"");
}
