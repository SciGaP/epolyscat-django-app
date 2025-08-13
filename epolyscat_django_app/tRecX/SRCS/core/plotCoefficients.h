// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PLOTCOEFFICIENTS_H
#define PLOTCOEFFICIENTS_H

#include <string>
#include <vector>
#include "tools.h"

class Coefficients;

/// \ingroup Plot

/// \brief base class for transforming and plotting Coefficients
class PlotCoefficients
{
protected:
    static double _tagIntervalDefault;
    std::string _acceptedTag;
    double _tagMin,_tagMax;
    double _tagInterval; /// minimal plot intervals (<0...not set)
    std::vector<std::string> _header;
    bool _overWrite;
public:
    PlotCoefficients():_acceptedTag("NONE"),_tagMin(-DBL_MAX),_tagMax(DBL_MAX),_tagInterval(_tagIntervalDefault),_overWrite(false){}
    /// transform to plot and write to file
    virtual void plot(const Coefficients & C, const std::string & File,
                      const std::vector<std::string> & Head=std::vector<std::string>(0),std::string Tag="",bool OverWrite=false) const=0;
    /// brief name of output file type
    virtual std::string briefName() const=0;

    ///< plot only for _tagMin <= tag _tagMax (if tag is convertible to double)
    double tagMin() const {return _tagMin;}
    double tagMax() const {return _tagMax;}

    static void setDefaultPlotInterval(double Interval) {_tagIntervalDefault=Interval;} ///< reset the default plot interval
    void setPlotInterval(double Interval, bool Override=false /** allow reseting previously set value */) const ;///< set interval for present plot
    bool acceptPlot(std::string Tag) const;

};

#endif // PLOTCOEFFICIENTS_H
