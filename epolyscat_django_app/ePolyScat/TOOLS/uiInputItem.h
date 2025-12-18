// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef UIINPUTITEM_H
#define UIINPUTITEM_H

#include <string>
#include <vector>

class ReadInput;

/// complete info on input as used in UI
class UiInputItem
{
    static void writeLines(std::fstream * Stream, const std::vector<UiInputItem> & Items);
protected:
    std::string _category;
    std::string _name;
    std::string _line;
    std::string _value;
    std::string _docu;
    std::string _defVal;
    std::string _flag;
    std::string _allow;
    std::string _texdocu;
    std::string _inpVal; // value as input (not superseeded by default)


public:
    static std::string markupGet(std::string Mark, std::string Line);    static std::string markUp(const std::string Mark,const std::string Value){return "\\"+Mark+"{"+Value+"}";}
    static void write(std::fstream * Stream, const std::vector<UiInputItem> & Items);

    UiInputItem(){};
    UiInputItem(std::string Line);
    UiInputItem(std::string Category,std::string Name,std::string Line,std::string Docu,std::string Default,std::string Flag, std::string Allow):
        _category(Category),_name(Name),_line(Line),_docu(Docu),_defVal(Default),_flag(Flag),_allow(Allow),_inpVal("")
    {if(_flag=="")_flag=Category+":"+Name;}
    std::string category() const{return _category;}
    std::string line() const{return _line;}
    std::string name() const{return _name;}
    std::string value() const{return _value;}
    std::string docu() const{return _docu;}
    std::string defaultValue() const{return _defVal;}
    std::string inputValue() const{return _inpVal;}
    std::string flag() const{return _flag;}
    std::string str() const {return _category+":"+_name;}
    int pos() const {return std::stoi(line());}

    std::string listDef() const; ///< string as item will appear in list-of-inputs
    std::string strMarkup() const; ///< full info in mark-up form (single line)
    /// adds LaTeX-like raw string (e.g. enter \latex{expression} as  R"mark(\latex{expression})mark")
    UiInputItem & texdocu(const char* RawText);
    bool matchesLinp(const std::string & LinpFile) const;

    UiInputItem & uiPage(std::string PageName) {return *this;}
    UiInputItem & inputUnits(std::string & UnitsString);

};

#endif // UIINPUTITEM_H
