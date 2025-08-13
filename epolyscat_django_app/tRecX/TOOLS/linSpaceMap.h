// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef MAPLINSPACE_H
#define MAPLINSPACE_H

#include <complex>
#include <vector>
#include "abort.h"

/// \ingroup Linalg
/// \brief map from a linear space onto itself (map is not necessarily linear)
///
/// V ...complies with LinSpaceVector
template <class V>
class LinSpaceMap{
    bool _allowAlias;
public:
    virtual ~LinSpaceMap(){}
    LinSpaceMap(bool AllowAlias=false):_allowAlias(AllowAlias){}
    /// Y <- A * Map(X) + B * Y
    virtual void apply(std::complex<double> A, const V & X,std::complex<double> B, V & Y) const=0;
    virtual bool applyAlias() const {return _allowAlias;} ///< true if &X==&Y is allowed in apply

    /// example object V (for use, e.g., in copy-constructors)
    virtual const V & lhsVector() const=0;
    virtual const V & rhsVector() const=0;

    /// update the map using a set of parameters
    virtual void update(double Time, const V* CurrentVector=0){};
};

#endif // MAPLINSPACE_H
