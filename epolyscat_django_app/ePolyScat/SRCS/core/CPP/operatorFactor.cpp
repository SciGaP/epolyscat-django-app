// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorFactor.h"

#include "readInput.h"
#include "tools.h"

using namespace std;

map<string,UseMatrix> OperatorFactor::matrix;

void OperatorFactor::readMatrix(ReadInput & Inp)
{
    string Name;
    unsigned int line=0;
    while (line<1000) {
        line++;
        Inp.read("FactorMatrix","name",Name,"","name to be used in operator definition",line);
        unsigned int nrow,ncol;
        Inp.read("FactorMatrix","nrow",nrow,"0","number of rows",line);
        Inp.read("FactorMatrix","ncol",ncol,tools::str(nrow),"number of columns (default == number of rows)",line);
        matrix[Name]=UseMatrix(nrow,ncol);
        string rowString;
        for(unsigned int i=0;i<nrow;i++){
            line++;
            vector<string> entry=tools::splitString(rowString,' ');
            if(entry.size()!=ncol)ABORT("too few entries for "+Name+", need row-wise BLANK SEPARATED list of all entries, is: "+rowString);
            for(unsigned int j=0;j<entry.size();j++)matrix[Name](i,j)=tools::string_to_double(entry[j]);
        }
        if(Inp.endCategory("FactorMatrix",line))return;
    }
    ABORT("FactorMatrix entry limited to below 1000 lines in total");
}
