// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef MAKEPLOT_H
#define MAKEPLOT_H

#include <string>

class ReadInput;
class Algebra;

class MakePlot
{
    std::string _file;
    unsigned int _nPoints;
    double _xMin,_xMax;
public:
    MakePlot(ReadInput & Inp);
    void addPlot(const Algebra & A) const;

    double xMin() const {return _xMin;}
    double xMax() const {return _xMax;}
    unsigned int nPoints() const {return _nPoints;}
    std::string file() const {return _file;}
};

#endif // MAKEPLOT_H
