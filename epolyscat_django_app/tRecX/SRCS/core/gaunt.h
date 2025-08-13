// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef GAUNT_H
#define GAUNT_H

#include "qtEigenDense.h"
#include "tools.h"

///@brief OBSOLESCENT (in new code, use Gaunts instead)
class Gaunt{
  int lmax;                                  // lmax of the calculation
  int order;                                 // 2*lmax+1
  std::vector<Eigen::MatrixXd> Passoc;       //Passoc[m](l,x)
  Eigen::VectorXd quadX,quadW;

public:
  static Gaunt main;

  Gaunt();
  void initialize(int lmax);
  double coeff(int lc,int l1,int l2, int mc,int m1,int m2);
  bool coeff_isZero(int lc,int l1,int l2, int mc,int m1,int m2);
  void test();

  double twoYintegrals(int l1,int l2,int m1,int m2, std::string operat);
};

#endif // GAUNT_H
