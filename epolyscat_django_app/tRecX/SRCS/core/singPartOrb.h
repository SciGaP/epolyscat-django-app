// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef SINGPARTFUNC_H
#define SINGPARTFUNC_H

#include <vector>
#include <string>

class UseMatrix;
class ReadInput;
class Gaussian;

class SingPartOrb
{
public:
    SingPartOrb();
    virtual ~SingPartOrb(){}

    /// values at a set of cartesian points
    virtual void val(std::vector<std::vector<double> > XYZ, UseMatrix & Val) const;
    virtual unsigned int size() const {return trans.cols();}
    virtual void readTrans(ReadInput & Inp, const std::string Type);

private:
    const Gaussian gauss;
    UseMatrix trans;
};

#endif // SINGPARTFUNC_H
