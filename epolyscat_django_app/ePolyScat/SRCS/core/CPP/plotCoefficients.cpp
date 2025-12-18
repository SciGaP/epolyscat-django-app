// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "plotCoefficients.h"

#include "printOutput.h"

// initialize to accepting any interval
double PlotCoefficients::_tagIntervalDefault=-1.;

void PlotCoefficients::setPlotInterval(double Interval, bool Override) const {
    if(Interval!=_tagInterval and _tagInterval>=0. and Override)
        PrintOutput::message(Str("reset plot interval value from")+_tagInterval+"to"+Interval);
    if(_tagInterval<0. or Override)const_cast<PlotCoefficients*>(this)->_tagInterval=Interval;
}

bool PlotCoefficients::acceptPlot(std::string Tag) const
{
    // if empty tag, plot is accepted unconditionally
    if(Tag=="")return true;

    // do not plot if tag did not change
    if(Tag==_acceptedTag)return false;

    //HACK: failure if conversion not covered
    if(not tools::doubleInside(tools::string_to_double(Tag),_tagMin,_tagMax))return false;

    // reject interval too short (unless no interval was set)
    if(_tagInterval>0.
            and _acceptedTag!="NONE"
            and tools::string_to_double(Tag)-tools::string_to_double(_acceptedTag)<_tagInterval*(1.-1.e-10)
            )return false;

    const_cast<PlotCoefficients*>(this)->_acceptedTag=Tag; // not so nice...
    return true;
}
