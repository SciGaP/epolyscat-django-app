// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "readInputRange.h"

#include <string>
#include "readInput.h"
#include "printOutput.h"
#include "tools.h"
#include "algebra.h"

using namespace std;

ReadInputRange::ReadInputRange(const ReadInput &Inp)
{
    for(int k=0;k<Inp.inputTable.size();k++){
        if(Inp.inputTable[k]._value.find("|")!=string::npos){
            _categ.push_back(Inp.inputTable[k]._category);
            _names.push_back(Inp.inputTable[k]._name);
            _line.push_back(Inp.inputTable[k].lineCat);

            int n;
            string inUnit=Inp.unitSystem,sLow,sUp;
            if(_names.back().find("(")!=string::npos)
                inUnit=tools::stringInBetween(_names.back(),"(",")");
            range(Inp.inputTable[k]._value,inUnit,sLow,sUp,n);

            double vlow=Algebra(sLow).val(0.).real();
            double vup=Algebra(sUp).val(0.).real();

            string outUnit=Inp.unitSystem;
            _lowVal.push_back(Units::convert(vlow,inUnit,outUnit));
            _upVal.push_back( Units::convert(vup, inUnit,outUnit));
            step.push_back((_upVal.back()-_lowVal.back())/n);
        }
    }
}

void ReadInputRange::print() const {
        PrintOutput::title("Parameter range");
        PrintOutput::paragraph();
        PrintOutput::newRow();
        PrintOutput::rowItem(" ");
        PrintOutput::rowItem("line");
        PrintOutput::rowItem("from");
        PrintOutput::rowItem("to");
        PrintOutput::rowItem("step");
        for(int k=0;k<_names.size();k++){
            string parUnit("DEFAULT_SYSTEM");
            if(_names[k].find("(")!=string::npos)parUnit=tools::stringInBetween(_names[k],"(",")");
            PrintOutput::newRow();
            PrintOutput::rowItem(_names[k]);
            PrintOutput::rowItem(_line[k]);
            PrintOutput::rowItem(Units::convert(_lowVal[k],"DEFAULT_SYSTEM",parUnit));
            PrintOutput::rowItem(Units::convert(_upVal[k],"DEFAULT_SYSTEM",parUnit));
            PrintOutput::rowItem(Units::convert(step[k],"DEFAULT_SYSTEM",parUnit));
        }
        PrintOutput::paragraph();
}

string ReadInputRange::low(const string Range){
    string units,lo,up;
    int n;
    range(Range,units,lo,up,n);
    if(units!="")lo+=" "+units;
    tools::cropString(lo);
    return lo;
}
void ReadInputRange::range(const string Range,string & Unit, string & Low,string & Up,int & N){
    Low=Range;Up=Range;N=0;
    if(Range.find("|")==string::npos)return;
    vector<string> valUnit=tools::splitString(Range,' ');
    if(valUnit.size()==0)return;
    vector<string> range=tools::splitString(valUnit[0],'|');
    Low=range[0];
    Up=range[1];
    N=1;
    if(range.size()==3)N=tools::string_to_int(range[2]);
    if(valUnit.size()>1)Unit=valUnit[1];
}
