// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEXDERIVED_H
#define INDEXDERIVED_H

#include "index.h"

/// \ingroup Index
/// \brief base class for derived Index's
class IndexDerived : public Index
{
public:
    IndexDerived();

    /// selected levels converted to grid
    IndexDerived(const Index * I, std::vector<unsigned int> Level, std::vector<unsigned int> Point=std::vector<unsigned int>(), std::vector<double> Limit=std::vector<double>());
protected:
    std::vector<const Index*> fromIndex; ///< list of indices from which the present is derived (if emtpy, identical Index tree structures assumed)
    std::string strNode() const;         ///< Index-specific data
};

#endif // INDEXDERIVED_H
