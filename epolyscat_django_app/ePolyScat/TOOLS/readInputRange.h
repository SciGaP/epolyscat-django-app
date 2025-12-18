// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef READINPUTRANGE_H
#define READINPUTRANGE_H

#include "multiParam.h"

class ReadInput;

/// \ingroup IO
/// @brief extract all parameter range inputs from Inp, arrange in MultiParam
///
/// for any real or integer parameter of given Category and Name,
/// an range can be specified in the format<br>
/// val0|val1|n<br>
/// all ranges will be put into a MultiParam that can be iterated
/// parameters read last will be steped through fastest
class ReadInputRange : public MultiParam
{
    std::vector<std::string> _categ,_names;
    std::vector<unsigned int> _line;
public:
    ~ReadInputRange(){}
    static std::string low(const std::string Range); ///< return "val0" of range "val0|val1|n"
    static void range(const std::string Range, std::string &Unit, std::string & Low, std::string & Up, int & N); /// decompose range string

    ReadInputRange(){}
    ReadInputRange(const ReadInput & Inp);
    void print() const;
    unsigned int line(unsigned int k) const {return _line[k];}
    std::string name(unsigned int k) const {return _names[k];}
};

#endif // READINPUTRANGE_H
