// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef READINPUTLIST_H
#define READINPUTLIST_H

#include <string>

#include "readInput.h"

/// \ingroup IO
/// \brief read input from ReadInput-generated "inputList"
class ReadInputList: public ReadInput
{
    std::string file;
public:
    static std::string itemDef(std::string Category, std::string Name, std::string Line){
        return Category+ReadInput::categoryTerminator+Name+"["+Line+"]=";
    }

    ReadInputList(std::string File);
    std::string readValue(const std::string Category, const std::string Name,
                          const std::string Default, const std::string Docu, unsigned int Line, std::string Flag, std::string Allow);
    /// return input that was actually input, do not substitute default
    std::string readInput(const std::string Category, const std::string Name, unsigned int Line);
    void write(const std::string Category, const std::string Name, unsigned int Line, std::string Value);
    void write(const std::string Flag, std::string Value);
};

#endif // READINPUTLIST_H
