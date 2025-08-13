// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include <string>
#include <map>
#include <vector>

#include "mpiWrapper.h"
#include "readInput.h"

size_t getPeakRSS();
size_t getCurrentRSS();
std::map<std::string, double> getMemInfo();

