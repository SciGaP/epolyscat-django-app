// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "discretizationCoulXNew.h"
#include "indexCoulX.h"

DiscretizationCoulXNew::DiscretizationCoulXNew(const Discretization *SurfDisc, double RC, double RMax, const std::vector<double> KGrid, bool BandOvr, bool PureBessel){
    name = "Rn.Phi.Eta.Rn";
    idx() = new IndexCoulX(SurfDisc->idx(),RC,RMax,KGrid,BandOvr,PureBessel);

//    hierarchy.push_back(Idx->hierarchy());
    this->continuityLevel.push_back(0);
    for(unsigned int i=0; i<SurfDisc->getAxis().size(); i++){
        axis.push_back(SurfDisc->getAxis()[i]); // possible not all needed axes?
    }
}
