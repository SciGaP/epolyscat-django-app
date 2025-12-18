// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef FILETOOLS_H
#define FILETOOLS_H

#include <vector>
#include<string>

namespace tools
{
std::string findFirstLine(std::string File, const std::vector<std::string> & String); ///< return the first line that contains any of the strings

}




#endif // FILETOOLS_H
