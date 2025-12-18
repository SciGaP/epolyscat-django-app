// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef INPUTITEM_H
#define INPUTITEM_H

#include <climits>
#include <string>

#include "uiInputItem.h"


class InputItem:public UiInputItem{
    friend class ReadInput;
    friend class ReadInputList;
    friend class ReadInputRange;

    InputItem(){}
    InputItem(std::string Category,std::string Name,unsigned int Line, std::string Flag,std::string Docu,std::string Default,std::string Allow):
        UiInputItem(Category,Name,std::to_string(Line),Docu,Default,Flag,Allow),lineCat(Line),machinePos(INT_MAX)
    {}

    unsigned int lineCat;
    unsigned int lineFile;
    unsigned int machinePos; // position in machine-generateable list (if any)
    bool admissible(); ///< check match of input value with values listed in docu as "some text: value1...explanation, value2...explantation"
public:
    std::string docuLine(int MinDots=20, int Indent=25, int LineBreak=60) const;
    InputItem &writeLaTeX(std::string Dir); ///< LateX-formated docu as category_name.tex-file in Dir
    InputItem(std::string Category,std::string Name)
        :UiInputItem(Category,Name,0,"-","dummy","",""),lineCat(0),machinePos(INT_MAX){}
};


#endif // INPUTITEM_H
