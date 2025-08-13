// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include <vector>
#include "basisMatNumbers.h"
#include "stringTools.h"
#include "readInput.h"
#include "algebra.h"
#include "units.h"

BasisMatNumbers::BasisMatNumbers(const Eigen::MatrixXcd Mat){_mat=Mat;}

static std::vector<std::complex<double> > line_to_vector(std::string Line){
    std::vector<std::string> parts=tools::splitString(Line,' ');
    std::vector<std::complex<double> >res;
    for(std::string s: parts)res.push_back(tools::string_to_complex(s));
    return res;
}

/// Matrix: name, kind, rowsXcols
///<br> ....
///<br> ....
///<br> MyFavoriteName, matrixElements, 5 x 3
///<br> a00 a01 a02
///<br> a10 a11 a12
///<br> ...
///<br> a40 a41 a42
///<br> MyDiagMatrix, diagonal, 1 2 3 4 5
///<br> AnotherName, anotherKind
///
/// Note: row-elements are blank-separated
BasisMatNumbers::BasisMatNumbers(ReadInput &Inp, int &Line){
    std::string shape;
    Inp.read("Matrix","shape",shape,"0x0","format \"rows x cols\" or \"d0 d1 d2 ...\"",Line);

    std::string units;
    Inp.read("Matrix","units",units,"noConversion","input units: noConversion, eV|energy",Line);
    std::vector<std::complex<double> > row;
    if(shape.find("x")!=std::string::npos){
        int rows,cols;
        rows=tools::string_to_int(tools::splitString(shape,'x')[0]);
        cols=tools::string_to_int(tools::splitString(shape,'x')[1]);
        _mat.resize(rows,cols);
        for (int Lin0=++Line;Line<Lin0+rows;Line++){
            row=line_to_vector(Inp.lineAt("Matrix",Line));
            if(int(row.size())!=cols)ABORT(Str("specify matrix as")+rows+"lines of blank-separated values\nfound"+row+"at line="+Line);
            if(units!="noConversion")row=Units::convert(row,units);
            _mat.row(Line-Lin0)=Eigen::Map<Eigen::VectorXcd>(row.data(),row.size());
        }
        Line--; // to last line read
    }
    else {
        row=line_to_vector(shape);
        if(row.size()==0)ABORT(Str("specify diagonal as blank-separated values\nfound \"")+shape+"\" at line="+Line);
        if(units!="noConversion")row=Units::convert(row,units);
        _mat=Eigen::MatrixXcd::Zero(row.size(),row.size());
        for(size_t k=0;k<row.size();k++)_mat(k,k)=row[k];
    }
}

