// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "fileTools.h"

#include <iostream>
#include <fstream>

namespace tools{

std::string findFirstLine(const std::string File, const std::vector<std::string> & String){
    std::ifstream f(File.c_str());
    std::string line;
    while(getline(f,line)){
        for(std::string s: String)
            if(line.find(s)!=std::string::npos)return line;
    }
    return "";
}
}
