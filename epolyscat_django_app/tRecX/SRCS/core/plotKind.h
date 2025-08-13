// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PLOTKIND_H
#define PLOTKIND_H

#include <vector>
#include <string>
#include <map>
#include <memory>


class ReadInput;

/// standard definitions for plots
class PlotKind{
    static std::map<std::string,std::map<std::string,std::shared_ptr<PlotKind>>> plotKinds;

protected:
    std::vector<std::string> _axes;
    std::vector<std::string> _use;
    std::vector<unsigned int> _points;
    std::vector<std::vector<double> > _bounds;
public:
    PlotKind(){}
    PlotKind(const std::vector<std::string> Axes,
             const std::vector<std::string> Use,
             const std::vector<unsigned int> Points,
             const std::vector<std::vector<double> > Bounds)
        :_axes(Axes),_use(Use),_points(Points),_bounds(Bounds){}

    PlotKind(const std::string Kind /** ReadInput category "Plot_"+Kind */,
             ReadInput & Inp);

    const std::vector<std::string> &axes() const{return _axes;}
    const std::vector<std::string> &use() const{return _use;}
    const std::vector<unsigned int> &points() const{return _points;}
    const std::vector<std::vector<double> > &bounds() const{return _bounds;}

    virtual std::vector<std::vector<double>> grid() const; ///< grid definitions for PlotKind
    virtual std::vector<std::vector<double>> weig() const; ///< weight definitions for PlotKind

    static const PlotKind *definePlot(const std::string Kind, const std::string Hierarchy);
    static std::string defaultKind(const std::string Hierarchy);

    unsigned int size(const std::string Axis) const;
    void resize(const std::string Axis, unsigned int NewSize);
};

/// Eta2-distribution grids equidistant in Theta1 (not Eta1) and energy (kRn*)^2/2 (not kRn*)
class PlotJAD: public PlotKind{
    std::vector<double> kGrid;
    std::vector<double> eta1Grid;
    std::vector<double> eta2Grid;
public:
    static std::string format();
    PlotJAD(std::string Def);
    std::vector<std::vector<double>> grid() const; ///< grid definitions for PlotJAD
    std::string strDef() const; ///< shorthand definition for file name

};

#endif // PLOTKINDS_H
