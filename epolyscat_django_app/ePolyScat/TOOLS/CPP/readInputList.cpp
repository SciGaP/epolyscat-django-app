// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "readInputList.h"
#include "mpiWrapper.h"
#include <stdio.h>      /* fopen, fputs, fclose, stderr */
#include <iostream>
#include <fstream>
#include "abort.h"
#include "stringTools.h"
#include "str.h"

ReadInputList::ReadInputList(std::string File)
    :file(File)
{
    std::ifstream stream(file.c_str());
    if(not stream.is_open())ABORT("could not open input file '"+file+"'");

    if(file.find("/")!=std::string::npos)_outDir=file.substr(0,file.rfind("/")+1);
    else                            _outDir="";

    std::string line;
    if(MPIwrapper::isMaster())
        while(getline(stream,line))allLines.push_back(line);
    MPIwrapper::Bcast(allLines,MPIwrapper::master());
    stream.close();
}

void ReadInputList::write(const std::string Category, const std::string Name, unsigned int Line, std::string Value){
    std::ofstream stream;
    stream.open(file.c_str(),std::ios_base::out|std::ios_base::app);
    stream<<itemDef(Category,Name,tools::str(Line))<<Value;
    stream.close();
}

void ReadInputList::write(const std::string Flag, std::string Value){
    std::ofstream stream;
    stream.open(file.c_str(),std::ios_base::out|std::ios_base::app);
    stream<<Flag+"="+Value;
    stream.close();
}

std::string ReadInputList::readValue(const std::string Category, const std::string Name, const std::string Default,
                                     const std::string Docu, unsigned int Line, std::string Flag, std::string Allow){
    std::string def=InputItem(Category,Name,Line,Flag,Docu,Default,Allow).listDef();
    for(std::string line: allLines){
        // machine-generated input is assumed to be legal - all checks are by-passed
         if(line.find(def)==0){
            line=line.substr(def.length(),line.find_first_of(" \\")-def.length());
            return tools::cropString(line);
        }
    }
    return ReadInput::notFound+"_"+Category+":"+Name;
}

std::string ReadInputList::readInput(const std::string Category, const std::string Name, unsigned int Line){
    std::string def=itemDef(Category,Name,tools::str(Line));
    for(std::string line: allLines){
        // machine-generated input is assumed to be legal - all checks are by-passed
         if(line.find(def)==0){
            std::string inp=UiInputItem::markupGet("inputValue",line);
            std::string dfl=UiInputItem::markupGet("default",line);
            std::string val=line.substr(def.length(),line.find_first_of(" \\")-def.length());
            val=tools::cropString(val);
            if(val!=dfl)return val;
            return inp!=dfl and inp!=ReadInput::notFound ? val : inp; // value is default and was not input from file
        }
    }
    return ReadInput::notFound+"_"+Category+":"+Name;
}
