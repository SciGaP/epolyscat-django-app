// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "debugInfo.h"

#include "coefficients.h"
#include "readInput.h"

static bool _destruct=true;
namespace debug_tools {
int rank=-1;

#ifdef _DEVELOP_
static int _lev=10;
#else
static int _lev=1;
#endif

bool destruct(){return _destruct;}
int verboseLevel(){return _lev;}
void setVerboseLevel(){
    ReadInput::main.read("Flag","verbose",_lev,tools::str(_lev),"0...minimal, 10...maximal",1,"verbose");
}

static std::string _section;
void setSection(std::string Section){_section=Section;};
std::string section(){return _section;}

void unsetDestruct(){_destruct=false;}

bool hasNan(const Coefficients *C){
    if(C->orderedData()){
        for(int k=0;k<C->size();k++)
            if(std::isnan(C->orderedData()[k].real()) or std::isnan(C->orderedData()[k].imag()))return true;
    }
    else
        for(int k=0;k<C->childSize();k++)
            if(hasNan(C->child(k)))return true;

    return false;
}
}
