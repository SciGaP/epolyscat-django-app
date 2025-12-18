// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "multiParam.h"
#include "readInput.h"
#include "printOutput.h"
#include "tools.h"

using namespace std;

MultiParam::MultiParam(const vector<double> &Up, const vector<double> &Low, const std::vector<double> &Step)
    :_lowVal(Low),_upVal(Up),step(Step)
{ construct();}

MultiParam::MultiParam(ReadInput &Inp, std::vector<string> &Names){read(Inp,Names);}

void MultiParam::read(ReadInput &Inp, std::vector<string> &Names){

    Names.clear();
    std::string name;
    double upIn,lowIn;
    int nSteps;
    unsigned int line=0;
    do {
        line++;
        Inp.read("ParameterRange","name",name,ReadInput::notFound,"name of the parameter",line);
        Inp.read("ParameterRange","max",upIn,"-infty","upper boundary",line);
        Inp.read("ParameterRange","nSteps",nSteps,"0","number of steps from lower to upper",line);

        double lowDef=0.;
        if(nSteps==0)lowDef=upIn;
        Inp.read("ParameterRange","min",lowIn,tools::str(lowDef,12),"lower boundary",line);

        if(name==ReadInput::notFound)break;

        if(name.find("(")!=string::npos){
            upIn= Units::convert( upIn,tools::stringInBetween(name,"(",")"));
            lowIn=Units::convert(lowIn,tools::stringInBetween(name,"(",")"));
        }

        Names.push_back(name);
        _lowVal.push_back(lowIn);
        _upVal.push_back(upIn);
        if(_lowVal.back()>_upVal.back())ABORT("negative parameter range min,max="+tools::str(lowIn)+","+tools::str(upIn));

        if(abs(lowIn-upIn)<1.e-12*max(1.,abs(lowIn)+abs(upIn)))step.push_back(1.);
        else step.push_back((_upVal.back()-_lowVal.back())/max(1,nSteps));
    } while (name!=ReadInput::notFound);
    _names=Names;
}

void MultiParam::print() const {
    PrintOutput::title("Parameter range");
    PrintOutput::paragraph();
    PrintOutput::newRow();
    PrintOutput::rowItem(" ");
    PrintOutput::rowItem("from");
    PrintOutput::rowItem("to");
    PrintOutput::rowItem("step");
    for(int k=0;k<_names.size();k++){
        PrintOutput::newRow();
        PrintOutput::rowItem(_names[k]);
        PrintOutput::rowItem(_lowVal[k]);
        PrintOutput::rowItem(_upVal[k]);
        PrintOutput::rowItem(step[k]);
    }
    PrintOutput::paragraph();
}


void MultiParam::construct()
{
    if(_lowVal.size()==0)
        for(unsigned int k=0;k<_upVal.size();k++)_lowVal.push_back(0.);
    if(_lowVal.size()!=_upVal.size())ABORT("sizes of lower and upper boundary vectors differ");
    if(step.size()==0)
        for(unsigned int k=0;k<_upVal.size();k++)step[k]=_upVal[k]-_lowVal[k];
    if(_lowVal.size()!=step.size())ABORT("sizes of step and boundary vectors differ");
}


/// increment multi-parameter, initialize=empty vector, final=return false and empty vector
/// NOTE: rightmost indices run fastest
bool MultiParam::next(std::vector<double> & Par){
    if(Par.size()==0){
        if(_lowVal.size()==0)return false;
        Par=_lowVal;
        return true;
    }
    for(int d=_lowVal.size()-1;d>=0;d--){
        if(Par[d]+step[d]<=_upVal[d]+step[d]*1.e-10){
            Par[d]+=step[d];
            for(unsigned int l=d+1;l<_lowVal.size();l++)Par[l]=_lowVal[l];
            return true;
        }
    }
    // if it cannot be incremented, reset and return false
    Par.clear();
    return false;
}
