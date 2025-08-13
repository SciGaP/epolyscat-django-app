// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "../latexDocu.h"

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>

#include "abort.h"
#include "tools.h"
#include "str.h"
#include "readInput.h"

static std::string defaultLaTeXdir="TEX/input/";
static std::map<std::string,std::string> _texdocu;
static std::map<std::string,std::string> _categorySort;


static std::string fileNameString(const std::string S){
    std::string s(S);
    for(char c: " @/()_")std::replace(s.begin(),s.end(),c,'_');
    return s;
}


void LatexDocu::categoryAdd(std::string Category, std::string Sort, std::string TexString,std::string Tutorials){
    std::string res=_categorySort[Category];
    if(res=="")_categorySort[Category]=Sort;
    std::vector<std::string> tuts=tools::splitString(Tutorials,',');
    std::string usage("--- no tutorial examples at present ---");
    if(tuts.size()){
        usage="{\\bf Usage examples: }";
        for(std::string t: tuts)
            usage+="\\nameref{docu:tutorial:"+t+"} ";
    }
    add("category_"+Category,
        "\\subsection*{\\textcolor{red}{\\lcode{"+Category+"}}}\n"
        +usage+"\\\\%texdocu\n"+TexString);
    if(not ReadInput::main.flag("generate-docu","re-write the LaTeX docu from code"))return;
    std::ofstream cat((defaultLaTeXdir+"/category_"+Category+".tex").c_str());
    cat<<"\\rule{0ex}{0ex}\\vspace*{3ex}\n";
    cat<<category(Category);
    cat<<"\\label{docu:category:"+Category+"}\n";
    cat<<get(Category);
    cat<<"% --- list of names ---\n";
}

std::string LatexDocu::category(std::string Category){
    std::string res="\\subsection*{\\textcolor{red}{\\lcode{"+Category+"}}}\n";
    if(_texdocu.count("category_"+Category))res=_texdocu["category_"+Category];
    return res;
}


void LatexDocu::add(std::string Label, std::string TexString){
    std::string res=_texdocu[Label];
    if(res=="")_texdocu[Label]=TexString;
    else if (res!=TexString)
        ABORT("inconsistent re-definition of 'Label'\n"
              " this: "+TexString
              +"\nprev: "+res);
}

std::string LatexDocu::get(std::string Label){
    return _texdocu.count(Label)?_texdocu[Label]:"";
}

static std::string truncateLine(std::string Line){
    size_t maxLen=70;
    std::string line(Line);
    for(char t: "\n, "){
        while(line.length()>maxLen){
            auto pos=line.rfind(t);
            if(pos==std::string::npos)break;
            line=line.substr(0,line.rfind(t));
        }
    }
    if(line.length()<Line.length())line+="\\ldots";
    return line.length()>maxLen-20?
                line : Line.substr(0,maxLen-3);

}

void LatexDocu::categoryTable(){
    // table consisting of 1st column category name, 2nd column, first box all entries, 2nd box: first line of tex-docu
    if(not ReadInput::main.flag("generate-docu","re-write the LaTeX docu from code"))return;
    std::map<std::string,std::string> sort=
    {
        {"Discretization","Axis,Absorption,BasisConstraint,PolarOffCenter,Polar2D,Constrain"}
        ,{"Operators","Operator,Dipole,Constant,Pot3d"}
        ,{"Eigenproblem","Eigen,Trace,FloquetAnalysis"}
        ,{"Time propagation","TimePropagation,Initial"}
        ,{"External Field","Laser,Pulse"}
        ,{"Spectra","Spectrum,SpectrumPlot,Surface,Source"}
        ,{"Output","Plot,PlotFunction,Output,Print,Title,CoefficientsWriter"}
        ,{"Controls","Flag,DEBUG,Checks"}
    };
    std::ifstream x(defaultLaTeXdir+"/table.tex");
    std::string line;
    std::map<std::string,std::string> categoryFile;

    while (std::getline(x,line)){
        if(line.find("% begin category")==0){
            std::string c;
            while(std::getline(x,line) and line.find("% end category")!=0){
                if(line.find("\\nameref{docu:category:")==0){
                    c=line.substr(23,line.find("}")-23);
                    categoryFile[c]=line+"\n";
                } else
                    categoryFile[c]+=line+"\n";
            }
        }
    }
    x.close();


    std::ofstream t(defaultLaTeXdir+"/table.tex");
    t<<"% --- will be updated automatically: only add/remove complete category entries ---\n";
    t<<R"tex(
       {\small
       \begin{longtable}{rl}
       )tex";
    t<<"\\hline\\hline\n";

    for(auto s: sort){
        t<<"{\\bf "<<s.first<<"}&\\\\ \\hline \n";
        auto cats=tools::splitString(s.second,',');
        for(std::string c: cats){
            t<<"% begin category\n";

            std::string entry="\\nameref{docu:category:"+c+"}&";
            if(_categorySort.count(c))entry+=truncateLine(_categorySort[c]);
            entry+="\\\\ \n";
            if(_texdocu.count("category_"+c)
                    and _texdocu["category_"+c].find("%texdocu\n")!=std::string::npos){
                // current entry
                t<<entry;
                std::string brief=_texdocu["category_"+c];
                brief=brief.substr(brief.find("%texdocu\n")+9);
                if(brief.length()>12)t<<"&"<<brief.substr(0,brief.find("\n",4))<<"\\\\ \n";
                t<<" \\hline\n";
            }
            else if(categoryFile.count(c)){
                // pre-exsiting entry
                t<<categoryFile[c];
            }
            else {
                // default entry
                t<<entry<<" \\hline\n";
            }
            t<<"% end category\n\n";
        }
    }
    t<<R"tex(\hline
       \end{longtable}
       })tex";
    t.close();
}

static bool inputItemLess(const InputItem & a, const InputItem & b){
    if (a.category() != b.category())return a.category() < b.category();
    return a.name() < b.name();
}

static bool inputCategoryLess(const std::pair<std::string,std::vector<InputItem>> & A, const std::pair<std::string,std::vector<InputItem>> & B){
    return A.first<B.first;
}

void LatexDocu::write(std::string LaTeXdir,std::vector<InputItem> inputTable){
    if(not ReadInput::main.flag("generate-docu","re-write the LaTeX docu from code"))return;

    ReadInput::main.texdocuCategoryAdd("Flag","generate-docu,recompute,showDiscretization,showMatrices,showOperators,"
                                              "printMatrices,noChecks,generate-linp,DEBUGfem, DEBUGold, DEBUGparallel",
                                       "Values that can ONLY by set by command line Flags","");
    if(LaTeXdir=="")LaTeXdir=defaultLaTeXdir;
    if(not folder::exists(LaTeXdir))
        folder::create(LaTeXdir);

    // select items with unique category and name
    std::map<std::string,std::vector<InputItem>> srtTable;
    std::string prevCat("NONE"),prevName;
    for (InputItem it: inputTable){
        if(it.defaultValue().find("OBSOLETE")<2)continue; // do not show obsolete inputs
        bool hasItem=false;
        for(auto nam: srtTable[it.category()])hasItem=hasItem or (nam.name()==it.name());
        if(not hasItem)srtTable[it.category()].push_back(it);
    }

    std::ofstream all((LaTeXdir+"/all.tex").c_str());
    all<<"% --- machine generated - do not edit ----\n";
    prevCat="NONE";

    std::fstream *cat(0);
    std::string startInput="\\input{input/";
    std::set<std::string> lines;
    for(auto &categ: srtTable){
        if (not _texdocu.count("category_"+categ.first))categoryAdd(categ.first,"","","");
        cat=new std::fstream((LaTeXdir+"/category_"+categ.first+".tex").c_str(),std::ios_base::app);
        std::vector<std::string> srt=tools::splitString(_categorySort[categ.first],',');
        for(auto &itm: categ.second){
            itm.writeLaTeX(LaTeXdir);
            std::string line=startInput+itm.category()+"_"+fileNameString(itm.name())+"}";
            *cat<<line<<"\n"<<std::flush;
        }
        delete cat;
    }
    LatexDocu::categoryTable();
}
