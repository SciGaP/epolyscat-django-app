// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DEBUGINFO_H
#define DEBUGINFO_H

#include <string>

class Coefficients;
namespace debug_tools{
bool destruct();
void unsetDestruct();
bool hasNan(const Coefficients* C);
int verboseLevel();//{return _lev;}
void setVerboseLevel();
void setSection(std::string Section);
std::string section();
extern int rank;
}
#endif // DEBUGINFO_H
