// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "uiInputItem.h"

#include "stdio.h"
#include <iostream>
#include <fstream>

#include "str.h"
#include "folder.h"
#include "readInputList.h"
#include "printOutput.h"

static std::string categorySeparator=":";
static std::string nameSeparator=",";



std::string UiInputItem::markupGet(std::string Mark, std::string Line) {
    std::string s("NO_"+Mark);
    size_t pos=Line.find("\\"+Mark);
    if(pos!=std::string::npos){
        s=Line.substr(pos+Mark.length()+1,Line.find("}",pos)-pos-Mark.length()-1);
        s=s.substr(s.find("{")+1);
    }
    return s;
}
std::string UiInputItem::listDef() const{
    return ReadInputList::itemDef(_category,_name,line());
}

std::string UiInputItem::strMarkup() const {
    std::string s=listDef();
    s+=_value;
    s+=" "+markUp("inputValue",_inpVal);
    s+=" "+markUp("default",_defVal);
    s+=" "+markUp("docu",_docu);
    s+=" "+markUp("allow",_allow);
     return s;
}

UiInputItem::UiInputItem(std::string Line)
{
    // old linp-style
    _category=Line.substr(0,Line.find(":"));
    _name=Line.substr(Line.find(":")+1,Line.find("[")-Line.find(":")-1);
    _line=Line.substr(Line.find("[")+1,Line.find("]")-Line.find("[")-1);
    _value=Line.substr(Line.find("=")+1,Line.find("\\")-Line.find("=")-1);

    if(_value.find("\\")!=std::string::npos){
        std::cout<<Line<<std::endl;
        std::cout<<"BAD: "<<_value<<" "<<Line.find("=")<<" "<<Line.find("\\")<<std::endl;
        exit(0);
    }

    // remove blanks
    while(_value[0]==' ')_value=_value.substr(1);
    while(_value.back()==' ')_value=_value.substr(0,_value.length()-1);


    // proper markup style
    _docu  =markupGet("docu",Line);
    _defVal=markupGet("default",Line);
    _allow =markupGet("allow",Line);

}

bool UiInputItem::matchesLinp(const std::string &LinpFile) const{
    if(not folder::exists(LinpFile))return true;
    return true;
    ReadInputList linp(LinpFile);
    std::string val=linp.readValue(category(),name(),defaultValue(),docu(),tools::string_to_int(line()),flag(),"");
    if(val!=value())
        PrintOutput::warning("["+category()+":"+name()+"] -"+flag()+"="+value()+" does not match ="+val+" on file "+LinpFile);
    return val==value();
}

void UiInputItem::write(std::fstream *Stream, const std::vector<UiInputItem> &Items){
    if(Items.size()==0)return;

    // group items in same category
    auto it=Items.begin();
    while(it!=Items.end() and it->category()==Items[0].category())it++;
    writeLines(Stream,std::vector<UiInputItem>(Items.begin(),it));

    // recursively write further groups
    write(Stream,std::vector<UiInputItem>(it,Items.end()));
}

/// write non-default Items according to UiLayout
void UiInputItem::writeLines(std::fstream *Stream, const std::vector<UiInputItem> &Items){
    if(Items.size()==0)return;
    auto itBeg=Items.begin();
    auto itEnd=Items.begin();
    while(itEnd!=Items.end()){
        while(itEnd!=Items.end() and itEnd->line()==itBeg->line())itEnd++;
        if(itBeg==Items.begin()){
            // write header
            *Stream<<std::endl<<itBeg->category()<<categorySeparator<<" "<<itBeg->name();
            for(auto it=itBeg+1;it!=itEnd;it++)*Stream<<itBeg->category()<<nameSeparator<<" "<<it->name();
        }

        // values in present line
        *Stream<<std::endl<<"     "<<itBeg->value();
        for(auto it=itBeg+1;it!=itEnd;it++){
            if(it->value()==it->defaultValue())*Stream<<itBeg->category()<<nameSeparator<<" "<<it->value();
            else                               *Stream<<itBeg->category()<<nameSeparator<<" ";
        }
    }
    *Stream<<std::endl;
}

UiInputItem &  UiInputItem::texdocu(const char* RawText){
    _texdocu=RawText;
    return *this;
}

UiInputItem &UiInputItem::inputUnits(std::string &UnitsString) {
    size_t pos=_value.find_first_of(" ~");
    UnitsString=(pos==std::string::npos)?"DEFAULT_SYSTEM":_value.substr(pos+1);
    return *this;
}









