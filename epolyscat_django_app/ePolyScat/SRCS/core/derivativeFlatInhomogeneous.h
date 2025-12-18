// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORINHOMOGENEOUSTSURFF_H
#define OPERATORINHOMOGENEOUSTSURFF_H

#include "tools.h"
//#include "operator.h"
#include "derivativeFlat.h"
#include "projectSubspace.h"

class TsurffSource;

/** \ingroup Structures */
/// provide the tSurff source term for partially ionized subregion

class DerivativeFlatInhomogeneous: public DerivativeFlat{
  std::vector<TsurffSource* > tS;
  const OperatorTree* _op; // for debug
public:
  DerivativeFlatInhomogeneous(const OperatorTree* Op, double ApplyThreshold, const DiscretizationSpectral *ProjectionDisc=0, std::vector<TsurffSource* > Source=std::vector<TsurffSource* >(0));
  DerivativeFlatInhomogeneous(const OperatorTree* Op, double ApplyThreshold,
                              std::shared_ptr<ProjectSubspace> Project, std::vector<TsurffSource* > Source=std::vector<TsurffSource* >(0));
  ~DerivativeFlatInhomogeneous();

  void update(double Time, const Coefficients* CurrentVec=0);
  void apply(std::complex<double> A, const Coefficients& X, std::complex<double> B, Coefficients &Y) const;
  void apply(std::complex<double> A, CoefficientsLocal *localX, std::complex<double> B, CoefficientsLocal &Y) const;
  const Index* idx() const;
};

#endif // OPERATORINHOMOGENEOUSTSURFF_H
