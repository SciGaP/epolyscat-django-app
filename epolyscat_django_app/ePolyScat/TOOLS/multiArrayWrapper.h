// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef MULTIARRAYWRAPPER_H
#define MULTIARRAYWRAPPER_H
#include <vector>
#include "qtEigenDense.h"
#ifdef _USE_BOOST_
#include <boost/multi_array.hpp>
typedef boost::multi_array<double,3> multiarray3d;
typedef boost::multi_array<Eigen::MatrixXcd,4> multiarray4EigenXcd;
typedef boost::multi_array<Eigen::MatrixXcd,3> multiarray3igenXcd;
#else
typedef std::vector<double> multiarray3d;
typedef std::vector<double> multiarray4EigenXcd;
typedef std::vector<double> multiarray3igenXcd;
#endif
#endif // MULTIARRAYWRAPPER_H
