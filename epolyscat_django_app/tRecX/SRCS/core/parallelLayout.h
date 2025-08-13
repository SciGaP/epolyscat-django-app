// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARALLELLAYOUT_H
#define PARALLELLAYOUT_H

#include <string>
#include <vector>
#include <map>

class Index;
class ReadInput;
class ParallelCross;

class ParallelLayout
{
    std::string _hier,_sort;
    std::vector<int> _permute;
    static bool smallerCross(const ParallelCross * A, const ParallelCross * B);
    std::map<std::string,unsigned int> _floorSetupHost; ///< setup host for floor
public:
    ParallelLayout(const Index* Idx, std::string Sorting);
    static void read(ReadInput & Inp);
    void sort(std::vector<ParallelCross*> &Cross); ///< sort that first index runs fastest (after permutation)
    void setFloorHosts(const Index* Idx);
    int floorHost(const Index* Floor);
};

#endif // PARALLELLAYOUT_H
