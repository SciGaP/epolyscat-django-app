// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef LATEXDOCU_H
#define LATEXDOCU_H

#include <string>
#include "readInput.h"

namespace LatexDocu {

/// generate LaTeX documentation from information in InputItems

/// create an overview "table.tex" with category, allowed list of names, and brief descriptoin
void categoryTable();

/// add long version of category description into file "category_"+Category"+.tex"
void categoryAdd(std::string Category, std::string Sort, std::string TexString,std::string Tutorials);

/// retrieve category description from internal storage
std::string category(std::string Category);

/// add docu item to internal storage, use Label (e.g. "category_Axis" or "Axis_name"), TexString must be unique
void add(std::string Label, std::string TexString);

/// get TexString stored under Label from internal storage ("" if Label not found)
std::string get(std::string Label);

/// add new items and extend category table with InputItems from current run
/// tex-files will be in LaTeXdir
void write(std::string LaTeXdir, std::vector<InputItem> inputTable);
};

#endif // LATEXDOCU_H
