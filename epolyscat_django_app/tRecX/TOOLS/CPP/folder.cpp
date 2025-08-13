// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
/// this is the Linux version ... and dirty
#include "folder.h"
#include <stdlib.h>     /* system, NULL, EXIT_FAILURE */
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;



#ifdef _WIN32
#include <boost/filesystem.hpp>
#endif



bool folder::create(const string & Name){
#ifdef _WIN32
	return boost::filesystem::create_directory(Name);
#else
    // Remark: This executes a clone syscall, which may take a while (memory
    // is copied) and is possibly interrupted by SIGPROF if compiled with -pg.
    // This creates an infinite loop. So... don't compile with -pg.
    //
	// success returns 0, but that is logical false in cpp: therefore "not"
	return !system(("mkdir " + Name).c_str());
#endif
}

bool folder::exists(const string & Name){
#ifdef _WIN32
	return boost::filesystem::exists(Name);
#else
    ifstream infile(Name.c_str());
    return infile.good();
#endif
}
