// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "qtEigenDense.h"
#include "index.h"
#include "basisDvr.h"

class TestIntegration
{
public:
    TestIntegration();
    Eigen::MatrixXd GetBasisCoeff(const Index* idx);
    double IntegrateMonomial(int power,double UpBound, double LowBound);
    double IntegratePolynomial(Eigen::VectorXd C, double UpBound, double LowBound, int powerofX); //int b x^power
    double IntegrateLowerTriangle(Eigen::VectorXd C1, Eigen::VectorXd C2,double UpBound,double LowBound);
    double IntegrateUpperTriangle(Eigen::VectorXd C1, Eigen::VectorXd C2,double UpBound,double LowBound);
    double IntegrateElementEE(const Index *idx, int i, const Index *jdx, int j, std::string type);
};
