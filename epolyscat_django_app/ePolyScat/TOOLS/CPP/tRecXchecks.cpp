// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "tRecXchecks.h"

#include <set>

#include "readInput.h"
#include "printOutput.h"
namespace tRecX{

std::vector<std::string> offList;

void read(){
    ReadInput::main.read("Checks","off",offList,"","do not perform checks as listed, if emty list switch of globally",1,"noCheck");
    if((offList.size()==0 and ReadInput::main.found("Checks","off")) or
            ReadInput::main.flag("noChecks","globally switch off all checks"))
            offList.insert(offList.begin(),1,"__GLOBAL__");
    if(offList.size())PrintOutput::DEVmessage(Sstr+"list of disabled checks: "+offList);
}
void print(){
    if(offList.size()==0)return;

    if(offList[0]=="__GLOBAL__"){
        PrintOutput::warning("checks switched off globally");
        return;
    }

    for(int k=0;k<offList.size();k++){
        PrintOutput::warning("check \""+offList[k]+"\" disabled");
    }
}

static std::set<std::string> doneList;
bool off(const std::string mess) {
    bool of;
    if(offList.size()==0)of=false;
    else if(offList[0]=="__GLOBAL__")of=true;
    else{
        of=false;
    for(int k=0;k<offList.size();k++)
        if(mess==offList[k]){
            of=true;
            break;
        }
    }
    if(of){
        doneList.insert(mess);
    }
    return of;
}

void setOff(std::string Off){
            offList.push_back(Off);
}

}
