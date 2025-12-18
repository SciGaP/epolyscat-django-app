// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisMatAbstract.h"

const UseMatrix BasisMatAbstract::useMat() const {
    return UseMatrix::UseMap(const_cast<BasisMatAbstract*>(this)->_mat.data(),_mat.rows(),_mat.cols());
}

const std::vector<const Eigen::MatrixXcd*> BasisMatAbstract::mats() const {
    return std::vector<const Eigen::MatrixXcd*>(1,&_mat);
}
