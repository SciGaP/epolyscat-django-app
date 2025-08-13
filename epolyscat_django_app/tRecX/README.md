# tRecX - time-dependent Recursive indeXing
tRecX = tSurff+irECS - a universal Schroedinger solver
Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)

This code is originally from the group of Armin Scrinzi at the LMU Munich,
see [here](https://trecx.physik.uni-muenchen.de/home.html).

Great thanks to Tobias Koelling, who did all the porting to Windows and Mac OS
and who initiated this documentation.

LICENSE:

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

# Publication

http://arxiv.org/abs/2101.08171

# Directories
 
    tRecX/SRCS        specific sources for the code
    tRecX/TOOLS       sources for wider use
    LIBRARIES         default libraries (replace by system libraries, where available)

# Compilation

1. Supply libraries, compile where needed (see details below)

2. Create your system specific Makefile, using the CMAKE utility

   prompt> cmake .

3. Compile by

   prompt> make


### Trouble shooting


1. prompt> make clean

2. Start over


### Needed libraries

The following libraries are needed:

* ARPACK (see comments below)
* LAPACK
* BLAS
* CBLAS
* FFTW
* boost-system

All libraries used are public domain and can be downloaded from the web.

Except for alglib, all other libraries are available in standard Linux distributions.

Warning: lapack is broken in some distributions, in that case compile manually.


### Required tools

* C compiler, e.g. gcc (reasonably new)
* C++ compiler, e.g. g++
* cmake
* make
* possibly git (for retrieval of the code from the git repository)


### Compiler specifics

* gcc 4.8.5 (e.g. shipped with SuSe leap): run with -std=gnu++11 (will not work with -std=c++11)

# Documentation

Overview of classes and code structure by Doxygen

## Build the Doxygen documentation

The documentation can be built using

    prompt> doxygen

This requires the following programs to be installed:

* ``doxygen`` (>= 1.8.3 works best) see [here](http://www.doxygen.org)
* ``dot`` for inheritance graphs see [here](http://www.graphviz.org/)

# Ports

The following ports may not be working any more (not maintained)

## Mac

The existing CMakeLists.txt works on some versions of Mac OS. If it does not,
you may try to edit the Mac-specific definitions in CMakeCache.txt
or examine the beginning of LIBRARY/CMake.include, mostly include and library paths,
for example:

if(CMAKE_SYSTEM_NAME MATCHES Darwin) <br/>
set (tRecX_LIBRARY_PATH /usr/local/Cellar/gcc/8.2.0/lib/gcc/8) <br/>
set (tRecX_INCLUDE_PATH /usr/local/Cellar/gcc/8.2.0/lib/gcc/8) <br/>
endif()

## Windows

### git

Install "Git for Windows" from [here](http://msysgit.github.io/). Maybe [Tortoise Git](https://code.google.com/p/tortoisegit/) is a helpful addition, it has a better context-menu integration for Windows Explorer.

### MinGW

It is possible to get all this running under Windows using [MinGW](http://www.mingw.org/wiki/Getting_Started). I did the following selection, however it may be possible to replace *Developer Toolkit* by *MSYS*:

* MinGW Compiler Suite
    * C Compiler
    * C++ Compiler
    * Fortran Compiler
* MinGW Developer Toolkit ?!

After the installation, ``wget`` has to be installed as well (on the command line, e.g. ``cmd``):

    mingw-get install msys-wget-bin

### MinGW-w64

There is a 64bit implementation of MinGW, which I expect to work as well. However, I did not do any tests. MinGW-w64 can be found [here](http://mingw-w64.sourceforge.net/).

### Compile

Do the preparations as mentioned above.

MinGW supplies two versions of make: ``mingw32-make`` in the MinGW package and ``make`` in MSYS. If you run from ``cmd``, use ``mingw32-make``, as ``make`` needs the MSYS-shell.


