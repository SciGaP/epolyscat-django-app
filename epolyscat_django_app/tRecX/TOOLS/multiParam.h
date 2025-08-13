// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef MULTIPARAM_H
#define MULTIPARAM_H
#include "toolsHeader.h"

class ReadInput;
///< a general multi-index class
class MultiParam{

protected:
    std::vector<std::string> _names;
    std::vector<double> _lowVal,_upVal,step;
    void construct();
    void read(ReadInput & Inp,std::vector<std::string> & Names);

public:
    MultiParam(){}

    /// \brief MultiParam parameter range defined on rectangle
    /// \param Up     upper bounderies
    /// \param Low    lower boundaries
    /// \param Step   step size for stepping through rectangle
    MultiParam(const std::vector<double> & Up, const std::vector<double> & Low=std::vector<double>(0), const std::vector<double> & Step=std::vector<double>(0));

    /// read multi-parameter ranges from input
    MultiParam(ReadInput & Inp,std::vector<std::string> & Names);


    /// increment multiple parameters, initialise=empty vector, final=return false and empty vector
    bool next(std::vector<double> & Par);
    std::vector<double> first(){return _lowVal;} //!< set first multi-index
    bool empty() const {return _lowVal.size()==0;}
    unsigned int size(){return _lowVal.size();}
    void print() const;
    double lowVal(int I){return _lowVal[I];}
    double upVal(int I){return _upVal[I];}
};

#endif // MULTIPARAM_H

