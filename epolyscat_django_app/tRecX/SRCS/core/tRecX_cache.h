// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef CHACHING_H
#define CHACHING_H

#include <string>

class ReadInput;


namespace tRecX_cache
{
void read(ReadInput & Inp); ///< set cache names and sources
std::string source(std::string Name); ///< source identifier (or "" if not in cache list)
std::string streamName(std::string Name); ///< stream name (or "" if not in cache list)
};

#endif // CHACHING_H
