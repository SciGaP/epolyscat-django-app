// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorZDxZD.h"
#include "qtEigenDense.h"
#include "tools.h"
#include "useMatrix.h"

using namespace std;
#include "eigenNames.h"

OperatorZDxZD::OperatorZDxZD(std::vector<const UseMatrix *> Mat, std::string Kind)
{
    if(Mat[0]->rows()!=Mat[0]->cols() or Mat[1]->rows()!=Mat[1]->cols())
        DEVABORT("only for square diagonal matrices: "+Kind);

    std::vector<std::complex<double>> d0,d1,dd;
    for(int k=0;k<Mat[0]->rows();k++)d0.push_back((*Mat[0])(k,k));
    for(int k=0;k<Mat[1]->rows();k++)d1.push_back((*Mat[1])(k,k));
    for(std::complex<double> v0: d0)
        for(std::complex<double> v1: d1)
            dd.push_back(v0*v1);
    construct(Eigen::Map<Eigen::VectorXcd>(dd.data(),dd.size()),Kind);
    facs.push_back(OperatorZD(Mat[0],Kind+"factor0"));
    facs.push_back(OperatorZD(Mat[1],Kind+"factor1"));
}
Eigen::MatrixXcd OperatorZDxZD::matrixFactor(int D) const{
    return facs[D].matrixFactor(0);
}
