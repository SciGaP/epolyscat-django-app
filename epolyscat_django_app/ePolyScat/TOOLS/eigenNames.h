// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef EIGENNAMES_H
#define EIGENNAMES_H

/// frequent eigen names (instead of "using namespace Eigen")
///
/// Eigen.3.3.4 contains Eigen::Index, which conflicts with tRecX's Index
using Eigen::Map;

// this may be not too good an idea
using Eigen::Matrix;
using Eigen::RowMajor;
using Eigen::Dynamic;

using Eigen::Vector2i;
using Eigen::Vector3i;
using Eigen::VectorXi;
using Eigen::RowVector3i;
using Eigen::RowVector4i;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVector3d;

using Eigen::MatrixXcd;
using Eigen::Vector3cd;
using Eigen::VectorXcd;
using Eigen::ArrayXcd;
using Eigen::RowVector3cd;
using Eigen::RowVectorXcd;


#endif // EIGENNAMES_H
