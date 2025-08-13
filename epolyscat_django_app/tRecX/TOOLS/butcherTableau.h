// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BUTCHERTABLEAU_H
#define BUTCHERTABLEAU_H

#include <map>
#include <vector>
#include <string>
#include <complex>
///@brief collection of various explicit Runge-Kutta methods

/// from: L. Lapidus and J. Seinfeld,
///  <br>     "Numericel Solution of Ordinary Differential Equations"
///  <br>     Academic Press, New York and London, 1971
class ButcherTableau{
    std::string _name;
    int _consistency;
    static std::map<std::string,ButcherTableau> tab;
    std::vector<std::vector<double> > _a;
    std::vector<double> _b,_c;
public:

    ButcherTableau():_consistency(0){}
    ButcherTableau(std::string Name);

    std::string name() const {return _name;}
    int consistency() const {return _consistency;}
    std::vector<std::vector<double> > a() const {return _a;}
    std::vector<double> b() const {return _b;}
    std::vector<double> c() const {return _c;}

};


#endif // BUTCHERTABLEAU_H
