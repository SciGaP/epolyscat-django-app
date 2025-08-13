// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISABSTRACT_H
#define BASISABSTRACT_H

#include <deque>
#include <string>
#include "labelled.h"

#include "basisFunction.h"

class ComplexScaling;
class Index;
class Axis;
class UseMatrix;
class BasisSetDef;
class Coordinate;

class BasisSub;
class BasisIntegrable;
class BasisGrid;
class BasisGridQuad;
class BasisNdim;
class BasisOrbital;
class BasisCI;
class BasisHybrid;

/** @defgroup Basissets Basis sets
 * \ingroup DiscretizationClasses
 * \brief One and multi-dimensional basis functions, grid; operations on basis
* @{
*/

/// @brief Abstract base for all basis sets and grids
class BasisAbstract: public Labelled<BasisAbstract>
{
protected:
    std::string _name;
public:
    virtual std::string name() const {return _name;}

    virtual unsigned int size() const=0;
    virtual std::string strDefinition() const; ///< complete definition of basis for use in factory

    virtual ~BasisAbstract(){}
    BasisAbstract(std::string Name="NONE"):_name(Name){}

    virtual bool operator==(const BasisAbstract &other) const;
    virtual std::string str(int Level=0) const; ///< print basis parameters
    virtual const BasisAbstract *remove(const std::vector<int> &RemoveK) const; ///< return BasisSub - basis set with subset of functions removed

    virtual bool isAbsorptive() const; ///< true if absorption acts on basis

    // override as needed...
    virtual bool isPeriodic() const {return false;}
    virtual bool isDVR() const {return false;}
    virtual bool isGrid() const{return false;}
    virtual bool isIndex() const{return false;}
    virtual bool isOrthonormal() const{return false;}

    static const BasisAbstract* factory(const BasisSetDef & Def);
    static const BasisAbstract* factory(const std::string & Def);
    static std::vector<const BasisAbstract*> select(std::function<bool(const BasisAbstract*)> Criterion);

    const BasisNdim* ndim() const;
    const BasisIntegrable* integrable() const;
    const BasisGrid* grid() const;
    const BasisGridQuad* gridQuad() const;
    const BasisSub* sub() const;
    const BasisOrbital* orbital() const;
    const BasisCI* ci() const;
    const BasisHybrid* hybrid() const;

    /// Return the physical quantity associated to the Index'th basis function
    /// This might be the magnetic/angular momentum, a grid point or sth similar.
    /// Meant to be used for plotting/visualization/purposes primarliy!
    virtual double physical(int Index) const{ return Index; }
};

/**
 *  @}
*/ // end Basissets

#endif // BASISABSTRACT_H
