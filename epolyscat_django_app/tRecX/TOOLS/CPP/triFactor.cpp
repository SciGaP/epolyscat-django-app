// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "triFactor.h"

///< detstroy the tri-diagonal factors (e.g. if parent matrix was changed)
void TriFactor::reset(){
    delete storeLU;
    storeLU=0;
    pivot.clear();
    delete shapeM;
    shapeM=0;
}


bool TriFactor::verify(const UseMatrix &M){
    UseMatrix vec=UseMatrix::Random(M.rows(),1);
    UseMatrix ve0;
    ve0=vec;
    solve('n',vec);
    ve0-=M*vec;
    if(ve0.maxAbsVal()>1.e12){
        ve0.transpose().print("error",0);
    }
    return ve0.maxAbsVal()<1.e12;
}
