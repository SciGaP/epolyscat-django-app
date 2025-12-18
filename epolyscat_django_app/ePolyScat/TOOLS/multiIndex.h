// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef MULTIINDEX_H
#define MULTIINDEX_H
#include "../TOOLS/toolsHeader.h"
#include "abort.h"

///< a general multi-index class
class MultiIndex{
public:
    MultiIndex(const std::vector<int> & M);
    /// increment multi-index, initial=empty vector, final=return false and empty vector
    /// NOTE: row-wise increment, i.e. rightmost indices run fastest
    bool next(std::vector<int> & i);
    bool nextCol(std::vector<int> & i); //!< column-wise increment, i.e. left-most indices run fastest
    void first(std::vector<int>& i); //!< set first multi-index
private:
    std::vector<int> m;
};

#endif // MULTIINDEX_H
