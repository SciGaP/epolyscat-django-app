// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef EIGENTOOLS_H
#define EIGENTOOLS_H

#include "qtEigenDense.h"
#include <string>


class UseMatrix;

/// \ingroup Linalg
/// \brief Convenience routines covering features missing in Eigen
namespace EigenTools
{
bool isIdentity(const Eigen::MatrixXcd & Mat, double Eps=0.);
bool isZero(const Eigen::MatrixXcd & Mat, double Eps=0.);
bool isSelfadjoint(const Eigen::MatrixXcd & Mat, double Eps=0.);
bool isSymmetric(const Eigen::MatrixXcd & Mat, double Eps=0.);
bool isAntisymmetric(const Eigen::MatrixXcd & Mat, double Eps=0.);
std::string str(const Eigen::MatrixXcd & Mat, int Digits=2, char Print='A', int Width=0);
Eigen::MatrixXcd &purgeInPlace(Eigen::MatrixXcd &Mat, double Eps=1.e-12);
Eigen::MatrixXcd purge(const Eigen::MatrixXcd &Mat, double Eps=1.e-12);
void saveAscii(std::string FileName,const UseMatrix &Mat);
void saveAscii(std::string FileName,const Eigen::MatrixXcd &Mat);
void readAscii(std::string FileName,Eigen::MatrixXcd & Mat);
bool compareMatrices(std::string FileA, std::string FileB, double Epsilon=1.e-12);
bool isNaN(const Eigen::MatrixXcd & M);
}

#endif // EIGENTOOLS_H
