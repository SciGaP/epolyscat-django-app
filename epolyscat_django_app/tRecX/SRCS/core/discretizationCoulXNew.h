// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONCOULXNEW_H
#define DISCRETIZATIONCOULXNEW_H

#include "discretizationDerived.h"

class DiscretizationCoulXNew:public DiscretizationDerived
{
public:
    DiscretizationCoulXNew(const Discretization *SurfDisc, double RC, double RMax, const std::vector<double> KGrid, bool BandOvr, bool PureBessel=false);
    OperatorAbstract*& getMapTo(std::string To){DEVABORT("cannot map to Parent");}
};

#endif // DISCRETIZATIONCOULXNEW_H
