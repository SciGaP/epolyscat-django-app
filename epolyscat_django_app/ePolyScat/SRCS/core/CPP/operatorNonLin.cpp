// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorNonLin.h"
//#include "operator.h"
#include "operatorSingle.h"
#include "operatorFloorNonLin.h"
#include "operatorMeanEE.h"

OperatorNonLin::OperatorNonLin(const std::string Name, const std::string &Definition, const Index *IIndex, const Index *JIndex)
    :OperatorTree(Name,Definition,IIndex,JIndex)
{
    childAdd(new OperatorTree(Name,Definition,IIndex,JIndex));
}
void OperatorNonLin::update(double Time, Coefficients* C){
    OperatorFloorNonLin* fNonLin;
    for(OperatorTree* op=firstLeaf();op!=nullptr;op=op->nextLeaf()){
        if((fNonLin=dynamic_cast<OperatorFloorNonLin*>(op->floor())))fNonLin->updateNonLin(Time,C);
    }
}
