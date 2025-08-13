// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISGRID_H
#define BASISGRID_H

#include "qtEigenDense.h"
#include <vector>
#include <complex>


#include "basisAbstract.h"

class BasisGrid: public BasisAbstract
{
    static std::map<std::vector<double>,const BasisGrid*>_allBasisGrid;
protected:
    BasisGrid(const BasisAbstract *Grid);
    BasisGrid(const std::vector<double> Mesh):BasisAbstract("grid"),_mesh(Mesh){}
    std::vector<double> _mesh;
public:

    static const BasisGrid* factory(const std::string Definition);
    static const BasisGrid* factory(const std::vector<double> Mesh);
    static const BasisGrid* factory(const BasisAbstract * Grid);

    bool isGrid() const{return true;}
    unsigned int size() const {return _mesh.size();}
    std::string strDefinition() const; ///< string fully defines Grid
    std::string str(int Level=0) const; ///< human-readable description of basis

    const std::vector<double> & mesh() const {return _mesh;}
    bool operator==(const BasisAbstract& Other) const;
    double upBound() const{return _mesh.back();}      ///< return upper boundary of interval
    double lowBound() const{return _mesh.front();}   ///< return lower boundary of interval
    double physical(int Index) const{return _mesh[Index];}
    Eigen::MatrixXcd mapInterpolate(const BasisAbstract* Bas, int Order=4) const;
};

#endif // BASISGRID_H
