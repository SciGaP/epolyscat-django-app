// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "plotSpectral.h"

#include "discretizationSpectral.h"
#include "index.h"
#include "operatorDiagonal.h"
#include "asciiFile.h"

using namespace std;

PlotSpectral::~PlotSpectral(){delete tmp;}

PlotSpectral::PlotSpectral(const DiscretizationSpectral * Spec):spec(Spec)
{
    tmp=new Coefficients(spec->mapFromParent()->iIndex);
}

void PlotSpectral::plot(const Coefficients &C, const std::string &File, const std::vector<std::string> &Head, string Tag, bool OverWrite) const{


    spec->mapFromParent()->apply(1.,C,0.,*tmp);

    std::vector<std::vector<double> > cols(4);
    intoCols(spec->eigenvalues(),tmp,cols);

    AsciiFile plotFile(File);
    vector<string>comm(Head);
    comm.push_back("");
    comm.push_back("Re(E)     Im(E)     Re(C)     Im(C)");
    plotFile.writeComments(comm);
    plotFile.writeCols(cols);
}

void PlotSpectral::intoCols(std::vector<std::complex<double>> Diag, Coefficients *C, std::vector<std::vector<double> > &Cols) const{
    if(C->isLeaf()){
        size_t pos=C->idx()->posIndex();
        for(unsigned int k=0;k<C->idx()->size();k++){
            Cols[0].push_back(Diag[k+pos].real());
            Cols[1].push_back(Diag[k+pos].imag());
            Cols[2].push_back(C->data()[k].real());
            Cols[3].push_back(C->data()[k].imag());
        }
    }
    else
        for(unsigned int k=0;k<C->childSize();k++)
            intoCols(Diag,C->child(k),Cols);
}
