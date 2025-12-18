// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
/** @defgroup DiscretizationClasses Discretization
 * \brief Kinds, coordinates, basis sets, indexing, coefficients
 */


/** @defgroup Discretizations Discretizaion kinds
 *  \ingroup DiscretizationClasses
 *  \brief Product of coordinate axes, spectral discretizations, grids etc.
 *  @{
*/



#ifndef __DISCRETIZATION__
#define __DISCRETIZATION__
#include "axis.h"
#include <memory>

// forward class delarations

class Wavefunction;
class Coefficients;
class Index;
class IndexFloor;
class IndexConstraint;
#include "axisTree.h"


/// \brief base class for all discretizations
class Discretization {
    Index* _idx;
public:
    friend class Coefficients;
    friend class OperatorData;
    friend class InFlux;
    friend class DiscretizationDerived;
    friend class DiscretizationConstrained;

    static bool indexSelect(const Index* I);

    Discretization():parent(0),_idx(0),constraint(0){}
    static Discretization * factory(ReadInput & Inp);

    /// Derived class destructor must delete Idx.
    virtual ~Discretization();

    /// alternative inverse overlap function for Coefficients
    virtual Coefficients& inverseOverlap(Coefficients& coeffs, Index* index, bool& CoeffsPerformInvOv);

    ///////////////////////////////////////////////////////////////
    // data
    std::string name;                   ///< name string (defaults to axis names if axes are given)
    std::vector<std::string> hierarchy; ///< strings naming the hierarchy (to be replaced by more structured info)
    Index*& idx(){return _idx;}              ///< the index system
    const Index * idx() const {return _idx;} ///< the index system
    std::vector<int> continuityLevel;   ///< (OBSOLESCENT - replace using Index::continuity(N)) hierarchy levels where continuity is imposed
    const Discretization * parent;      ///< derived discretization has a parent, top parent==0
    const IndexConstraint* constraint;  ///< Constraints imposed on the discretization
    ///////////////////////////////////////////////////////////////


    void show() const;                       ///< overview of the discretization
    std::string str(unsigned int Brief=0) const;                       ///< overview of the discretization
    void print(std::string File="", std::string Title="") const;   ///< print to file (default to cout)
    const std::vector<Axis> & getAxis() const { return axis;} ///< return the Discretization's axis
    const AxisTree* axisTree() const { return _axisTree.get(); } ///< return the Discretization's axis


    virtual std::string coordinates() const;

protected:
    std::vector<Axis> axis;
    std::shared_ptr<AxisTree> _axisTree;
    void setAxis(ReadInput &In,std::string Subset=""); //!< completely defines a discretization
    void construct();            //!< construct once axes are set up

//    const BasisSet* blockBasis(unsigned int LAxis, unsigned int N, const std::vector<const Index *> &Path) const; //!< return updated basis definition for current block

};
/** @}/ */

#endif
