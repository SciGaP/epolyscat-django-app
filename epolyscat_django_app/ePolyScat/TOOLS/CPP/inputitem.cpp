// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "inputitem.h"

#include <fstream>
#include <sstream>

#include "tools.h"
#include "readInput.h"

static std::string fileNameString(const std::string S){
    std::string s(S);
    for(char c: " @/()_")std::replace(s.begin(),s.end(),c,'_');
    return s;
}

static std::string linebreak(const std::string Line, int MaxLength, std::string Indent){
    std::string l(Line);
    std::replace(l.begin(),l.end(),'\n',' ');
    std::vector<std::string> parts=tools::splitString(l,' ');
    int n=0;
    std::string s;
    for(std::string p: parts){
        if(n>0 and int(n+p.length())>MaxLength){
            s+="\n"+Indent;
            n=Indent.length();}
        s+=p+" ";
        n+=p.length()+1;
    }
    return s;
}

InputItem& InputItem::writeLaTeX(std::string Dir){
    std::ofstream doc((Dir+"/"+category()+"_"+fileNameString(name())+".tex").c_str());
    doc<<"% --- machine generated - do not edit ----\n";
    doc<<"\\subsubsection*{\\textcolor{blue}{\\lcode{"+name()+"}}}\n";
    doc<<"\\label{docu:"+category()+ReadInput::categoryTerminator+fileNameString(name())+"}\n";
    doc<<"\\begin{lstlisting}\n";
    doc<<"  "<<linebreak(docuLine(5,15),75,"     ")<<"\n";
    doc<<"\\end{lstlisting}\n";
    doc<<"\\noindent"+_texdocu;
    doc<<"\\rule{0cm}{0cm}\\\\{\\small Go to \\nameref{sec:overview}}\\vspace*{2ex}\n";
    return *this;
}

std::string InputItem::docuLine(int MinDots, int Indent, int LineBreak) const{
    std::stringstream* doc=new std::stringstream();
    *doc << _name;
    std::string defFlag=_category+":"+_name;
    if(flag()!=defFlag)*doc<<" (-"+flag()+")";
    std::string defVal=_defVal;
    if(defVal==tools::str( DBL_MAX))defVal="+inf";
    if(defVal==tools::str(-DBL_MAX))defVal="-inf";
    *doc<<" ["+defVal+"]..";
    // fill up with dots...
    int d=_name.length()+defVal.length();
    if(_flag!=defFlag)d+=_flag.length()+4;
    for(;d<MinDots;d++)*doc<<".";

    std::string docu=_docu;
    size_t allowAll=docu.find("{allow all}");
    if(allowAll!=std::string::npos)docu.erase(docu.begin()+allowAll,docu.begin()+allowAll+std::string("{allow all}").length());
    std::vector<std::string>allow,dum;
    if(docu.find(":")!=std::string::npos and std::count(docu.begin(),docu.end(),',')>2){
        tools::splitString(docu.substr(docu.find(":")+1),",",allow,dum,"[({","])}");
        docu=docu.substr(0,docu.find(":")+1);
    }
    *doc<< docu;
    // line-break long docu's
    int len=0;
    std::string sp=" ";
    for(size_t k=0;k<allow.size();k++){
        len+=allow[k].length();
        if(len>LineBreak){
            *doc<<",\n"<<std::string(Indent,' ');
            len=0;
            sp="";
        }
        *doc<<sp<<tools::cropString(allow[k]);
        sp=",  ";
    }
    return doc->str();
}

bool InputItem::admissible(){
    std::string s(_docu);
    if(_docu.find("{allow all}")!=std::string::npos)return true;
    size_t cpos=tools::findFirstOutsideBrackets(_docu,":",ReadInput::leftBracket,ReadInput::rightBracket);
    if(cpos==std::string::npos)return true; // no specification, anything goes
    if(_value==ReadInput::noItem)return true; // no checking for dummy inputs
    if(_value==_defVal)return true; // default is always admissible
    std::string v=_value.substr(0,_value.find_first_of("[0123456789"));
    size_t vpos=s.find(v);
    if(vpos==std::string::npos)return false;
    if(vpos<cpos)return false;
    if(vpos-1!=cpos and s[vpos-1]!=' ' and s[vpos-1]!=',')return false;
    if(tools::cropString(s.substr(vpos,s.find_first_of("[.,0123456789",vpos)-vpos))!=v) return false;
    return true;
}
