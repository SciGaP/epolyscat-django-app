// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "vectorValuedFunction.h"
#include "tools.h"
#include "printOutput.h"

std::map<std::string,std::shared_ptr<const VectorValuedFunction> > VectorValuedFunction::_list;
bool VectorValuedFunction::validArguments(std::vector<double> X, std::string Message) const{

    double eps=1.e-12;
    bool valid=true;
    for(size_t k=0;k<_lowLim.size();k++)valid=valid and X[k]<_lowLim[k]-eps;
    for(size_t k=0;k<_upLim.size();k++) valid=valid and X[k]>_upLim[k]+eps;
    if(not valid and Message!="")abortInvalidArguments(tools::str(X)+" "+Message);
    return valid;
}

std::string VectorValuedFunction::info(int Number) const {
    std::string s=name()+" ("+coordinates()+") ["+tools::str(length())+"]";
    if(Number>-1)s=s+" No."+tools::str(Number);
    return s;
}

std::string VectorValuedFunction::list(){
    return tools::listMapKeys(_list,", ");}


void VectorValuedFunction::print(){
    if(_list.size()==0)return;

    PrintOutput::title("Vector-valued functions by name");
    for(auto f: _list){
        PrintOutput::newLine();
        PrintOutput::lineItem(f.first,f.second->info());
    }
}
