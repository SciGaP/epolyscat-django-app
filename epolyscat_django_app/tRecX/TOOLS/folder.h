// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef FOLDER_H
#define FOLDER_H

#include <string>

/// \ingroup IO
/// \brief check files and folders
namespace folder
{
bool create(const std::string & name);
bool exists(const std::string & name);
}

#endif // FOLDER_H
