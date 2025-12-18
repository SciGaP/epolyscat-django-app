// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef COORSYSTEM_H
#define COORSYSTEM_H

#include <algorithm>
#include <vector>
#include <string>
#include "qtEigenDense.h"

/// \ingroup Coordinates
/// \brief abstract base class for coordinate systems
class CoorSystem
{
    std::vector<std::string> standardCoors(std::string Name) const; ///< return standard coordinate names and sorting
    static CoorSystem* trySystem(CoorSystem* &S);

protected:
    std::string _standardName; ///< standard coordinate names and sorting, e.g. Phi.Eta.R for 3d polar
    std::string _name; ///< actual name and sorting, e.g. Phi1.R1.Eta1
    std::string _ref;  ///< reference coordinate system (e.g. X.Y.Z)

    CoorSystem(std::string StandardName, std::string RefCoor,std::string Name="");

    virtual std::vector<double> _toRef(const std::vector<double> & Coor) const =0;      ///< to reference coordinate system
    virtual std::vector<double> _fromRef(const std::vector<double> & Ref) const =0;     ///< from reference coordinate system
    virtual std::vector<double> _jacRefdCoor(const std::vector<double> & Coor) const =0;///< dim x dim jacobian matrix: d(Ref)/d(Coor) at Coor

public:
    virtual ~CoorSystem(){}
    static CoorSystem * factory(std::string System);
    static CoorSystem * factory(std::string OutSys,std::string InSys);

    int dim() const { return std::count(_name.begin(),_name.end(),'.')+1;} ///< dimension = number of single coordinates
    std::string name() const {return _name;} ///< coordinate system as input
    std::string refSystem() const {return _ref;} ///< reference coordinate system

    virtual std::vector<double> toRef(const std::vector<double> & Coor) const;      ///< to reference coordinate system
    virtual std::vector<double> fromRef(const std::vector<double> & Ref) const;     ///< from reference coordinate system

    /// Jacobian matrix: d(Ref)/d(Coor) evaluated at Coor
    ///
    /// ( dRef1/dCoor1   dRef1/dCoor2   dRef1/dCoor3 ...)<br>
    /// ( dRef2/dCoor1   dRef2/dCoor2   dRef2/dCoor3 ...)<br>
    /// ( dRef3/dCoor1   dRef3/dCoor2   dRef3/dCoor3 ...)<br>
    /// (                 ...                       )
    ///
    /// i.e. i'th row is dRef_i/dCoor_j, j=1,2,...
    ///
    /// (NOTE: CoordinateTrans returns transposed of this)
    Eigen::MatrixXd jacRefdCoor(const std::vector<double> & Coor) const;

    std::vector<int> posInStandard() const; ///< positions of coordinates in standard, e.g. StandardName=Phi.Eta.R -> Name=Phi.R.Eta gives {0,2,1}
};


#endif // COORSYSTEM_H
