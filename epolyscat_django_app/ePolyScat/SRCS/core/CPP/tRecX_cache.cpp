// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "tRecX_cache.h"


#include <map>
#include "readInput.h"
#include "folder.h"
#include "tools.h"
#include "str.h"

#include <istream>

namespace tRecX_cache
{
std::map<std::string,std::string> sources;
std::map<std::string,std::string> _alias={
    {"H0","hamiltonian"},
};

static std::string _dir="CACHE_TRECX/";

static std::string alias(std::string Name){
    return _alias.count(Name)?_alias[Name]:Name;
}

std::string source(std::string Name){
    if(not folder::exists(_dir))folder::create(_dir);
    return sources.count(alias(Name))?(_dir+sources[alias(Name)]+"_"+alias(Name)):"";
}

void read(ReadInput & Inp){
    std::string nameSource;
    Inp.read("Cache",ReadInput::anyName,nameSource,"",
             "blank-separated list --"
             " name1[source1] name2{[source2]} etc: "
             "hamiltonian...field-free Hamiltonian, source1 will be used as default for further sourceX",1,"cache")
            .texdocu(R"tex(
                     Objects identified by nameX can be stored and later retieved from disk cache.
                     \lcode{sourceX} is a freely chosen name for saving and retrieving the object identified by \lcode{name}.
                     Some consistency checks for the retrieved object will be made.
                     )tex");

    if(nameSource=="")return;

    std::vector<std::string> nams=tools::splitString(nameSource,' ');
    std::string srce("");
    for(std::string nam: nams){
        nam=tools::cropString(nam);
        if(nam.length()==0)continue;

        size_t beg=nam.find("["),end=nam.find("]");
        if(srce==""){
            if(beg==std::string::npos or end==std::string::npos)
                ABORT("must specify source for first name, have: "+nameSource);
            srce=nam.substr(beg,end-beg);
        }
        std::string anam=alias(nam.substr(0,beg));
        if(beg!=std::string::npos)sources[anam]=nam.substr(beg+1,end-beg-1);
        else                      sources[anam]=srce;
    }
};

}
