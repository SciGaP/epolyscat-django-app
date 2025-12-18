// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "readInput.h"
#include <stdio.h>      /* fopen, fputs, fclose, stderr */
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include "float.h"
#include  <algorithm>
#include "tools.h"
#include "folder.h"
#include "printOutput.h"
#include "algebra.h"
#include "mpiWrapper.h"
#include "readInputList.h"
#include "readInputRange.h"
#include "latexDocu.h"
#include "runDir.h"
#include "reuseRun.h"

using namespace std;
using namespace tools;

// digits in run numbering
int digitsRun(4);
std::string defaultOutputLevel;

static std::string defaultLaTeXdir="TEX/input/";

// special input values
const string ReadInput::noDefault = "NO_VALUE";
const string ReadInput::notFound = "NOT_FOUND";
const string ReadInput::empty = "READINPUT_EMPTY";
const string ReadInput::flagFound = "FLAG_FOUND";
const string ReadInput::flagOnly="FLAG_ONLY";
const string ReadInput::noItem = "ReadInput-noItem";
const string ReadInput::doNotAddToItemTable= "ReadInput-doNotAddToItemTable";
const string ReadInput::anyName = "ReadInput-anyName";

// file names and extensions
const string ReadInput::inputExtension = ".inp";
const string ReadInput::inputCopy = "inpc";
const string ReadInput::inputList = "linp";
const string ReadInput::docExtension = ".doc";

// syntax characters
const string ReadInput::itemSeparator = ",";
const string ReadInput::categoryTerminator = ":";
const string ReadInput::whiteSpace =" ";
const string ReadInput::comments = "!#*/";
const string ReadInput::quote = "'";
const string ReadInput::forbidden="`‘";
const string ReadInput::macro="#define";
const string  ReadInput::leftBracket=ReadInput::quote+"([{<|";
const string ReadInput::rightBracket=ReadInput::quote+")]}>>";
const string nonAlphaNumeric=ReadInput::leftBracket+ReadInput::rightBracket+ReadInput::itemSeparator+ReadInput::comments;
static string notInCategory=nonAlphaNumeric;

unsigned int ReadInput::mpiMaster=0;

static std::string markerSummary="Summary of input categories";

ReadInput ReadInput::main;

bool ReadInput::noInputFile() const {return _root=="";}

/// open and read the input file
ReadInput::ReadInput(int NCom, char *Com[]):_writeLinp(true){
    if(NCom<1)ABORT("need input file as first command line argument");
    string file=Com[1];
    //    if(file.length()<inputExtension.length() or
    //            file.rfind(inputExtension)!=file.length()-inputExtension.length())file+=inputExtension;
    construct(file,NCom,Com,true);
}

ReadInput::ReadInput(const string &File, int NCom, char *Com[], bool FromCommand):unitSystem("UNDEFINED"),_writeLinp(true)
{construct(File,NCom,Com,FromCommand);}

double ReadInput::readDouble(string Category, string Name, string Default, string Docu, unsigned int Line, string Flag, string Allow){
    double res;
    read(Category,Name,res,Default,Docu,Line,Flag,Allow);
    return res;
}
bool ReadInput::readBool(string Category, string Name, string Default, string Docu, unsigned int Line, string Flag, string Allow){
    bool res;
    read(Category,Name,res,Default,Docu,Line,Flag);
    return res;
}
int ReadInput::readInt(std::string Category,std::string Name,std::string Default,std::string Docu,unsigned int Line,std::string Flag,std::string Allow){
    int res;
    read(Category,Name,res,Default,Docu,Line,Flag,Allow);
    return res;
}

std::string ReadInput::nextRunDir(std::string Root){
    RunDir d(Root,digitsRun);
    if(d.nextUnused().dir()=="")ABORT("no more run numbers available in "+Root+"/");
    return d.dir();
}

void ReadInput::rootAndDir(std::string File){
    if(File.find(inputExtension)!=string::npos and File.rfind(inputExtension)+inputExtension.length()==File.length()){
        _outDir="NONE";
        _root=File.substr(0,File.rfind(inputExtension));
    }
    else if(File.find(inputCopy)!=string::npos and File.rfind("/"+inputCopy)+inputCopy.length()+1==File.length()){
        _outDir=File.substr(0,File.rfind(inputCopy));
        _root=_outDir.substr(0,_root.rfind(inputCopy));
        while(_root.back()=='/')_root.pop_back();
        _root=_root.substr(0,_root.rfind("/"));
    }
    else
        ABORT("input file name format must be either name"+inputExtension+" or root/0123/"+inputCopy);
}

std::string ReadInput::getInputFile(const std::string &File, int NCom, char **Com, bool FromCommand){
    string file(File);
    if(File!=flagOnly and FromCommand and NCom>1 and Com[1][0]!='-')
        if(File=="")file=Com[1]; // if second string on command line is not a flag, interprete as input file

    //HACK for AMP Gateway: assume .inp input file, but replace with restart run number
    for(int k=0;k<NCom;k++){
        if(std::string(Com[k]).find("-AMPrestart")==0){
            if(file.rfind(".inp")==file.length()-4)file=file.substr(0,file.length()-4);
            std::string from=Com[k];
            from=from.substr(from.rfind("=")+1);
            file=file+"/"+from;
            if(not folder::exists(file))
                ABORT("Restart: no input file \""+file+"\"\nUsage>tRecX InputName -restart=0123 for restart from InputName/0123");
        }
    }

    // construct input file name (presently inputCopy="inpc", inputExtension=".inp")
    std::string triedFile=file+" (input)";
    if(not ((file.find(inputCopy)!=string::npos  and file.rfind("/"+inputCopy)+inputCopy.length()+1==file.length()) or
            ((file.find(inputExtension)!=string::npos and file.rfind(inputExtension)+inputExtension.length()==file.length()))))
    {
        if(file[file.length()-1]=='/')file=file.substr(0,file.length()-1);
        triedFile+=", "+file+"/"+inputCopy+", "+file+inputExtension;
        if(file.rfind("/")+digitsRun+1==file.length()
                and file.find_first_not_of("0123456789",file.rfind("/")+1)==string::npos
                and folder::exists(file+"/"+inputCopy))file+="/"+inputCopy;  // is name/0123/inpc
        else if(folder::exists(file+inputExtension))file+=inputExtension; // is name.inp
    }

    if(file=="" and NCom==1){
        std::vector<std::string> commLin;
        for(int k=0;k<NCom;k++)commLin.push_back(std::string(Com[k]));
        showHelp(commLin);
    }
    return file;
}

void ReadInput::construct(const string & File, int NCom, char * Com[], bool FromCommand){

    // format category[specification]: is allowed
    notInCategory=nonAlphaNumeric;
    size_t pos;
    while((pos=notInCategory.find(']'))!=std::string::npos)notInCategory.erase(notInCategory.begin()+pos);
    while((pos=notInCategory.find('['))!=std::string::npos)notInCategory.erase(notInCategory.begin()+pos);


    for (int k=0;k<NCom;k++)commandLine.push_back(Com[k]);

    std::string file(getInputFile(File,NCom,Com,FromCommand));

    // construct name for input documentation file
    mainProgram=file.substr(0,file.find(inputExtension));
    if(NCom>0)mainProgram=Com[0];

    if(file=="" or file==inputExtension)showHelp(commandLine);

    if(file!=inputExtension and file!=flagOnly){
        ifstream stream(file.c_str(),std::ios::in);
        if(not stream.is_open()){
            PrintOutput::warning("EXIT: could not open input file '"+file+"'");
#ifdef _DEVELOP_
            std::abort();
#else
            std::exit(0);
#endif
        }

        string line;
        allLines.clear();
        if(MPIwrapper::isMaster())
            while(getline(stream,line))allLines.push_back(line);

        MPIwrapper::Bcast(allLines,MPIwrapper::master());
    }
    _outPrefix="";

    for(int n=0;n<NCom;n++)com.push_back(Com[n]);
    if(NCom==0){
        cout<<"--- flag input not enabled with "+file+" ---"<<endl;
    }

    rootAndDir(file);

    // get input macros and replace in allLines
    inputMacros();

    read("RunDirectory","digits",digitsRun,"4","how many digits for run numbering, min=4");
    digitsRun=std::max(digitsRun,4);

    std::string level;
    read("RunDirectory","outputLevel",defaultOutputLevel,"full","output level: full, low, off");
    PrintOutput::outputLevel(defaultOutputLevel);
}


std::string removeQuotesAndBlanks(string & ReturnString){
    // remove leading and trailing white spaces, outermost quotes
    ReturnString=tools::cropString(ReturnString);
    ReturnString=tools::unquote(ReturnString);
    ReturnString=tools::cropString(ReturnString);
    return ReturnString;
}


std::string ReadInput::strAdmissibleInputs(std::string Docu) const {
    if(Docu.find(":")==string::npos or Docu.find("{allow all}")!=string::npos)return " --- no input restrictions ---";
    Docu=Docu.substr(Docu.find(":")+1);
    while(tools::findFirstOutsideBrackets(Docu,",","[(","])")!=string::npos)
        Docu.replace(tools::findFirstOutsideBrackets(Docu,",","[(","])"),1,"\n");
    return Docu;
}

InputItem & ReadInput::findMatching(InputItem & Cur) {
    for(size_t k=0;k<inputTable.size();k++){
        InputItem & ii=inputTable[k];
        if(ii._category==Cur._category and ii._name==Cur._name and ii.lineCat==Cur.lineCat and ii._value==Cur._value){
            return ii;
        }
    }
    return Cur;
}

size_t findInput(const std::vector<std::string> & Lines, std::string Category, std::string Name){
    for (const std::string & l: Lines){
        if(l.find(Category+ReadInput::categoryTerminator)==0){
            if(l.find(Name,l.find(ReadInput::categoryTerminator))!=std::string::npos)return &l-&Lines[0];
        }
    }
    return string::npos;
}
bool duplicateInput(const std::vector<std::string> & Lines, std::string Category, std::string Name){
    size_t first=findInput(Lines,Category,Name);
    if (first==string::npos)return false;
    size_t second=findInput(std::vector<std::string>(Lines.begin()+first+1,Lines.end()),Category,Name);
    if(second==string::npos)return false;
    PrintOutput::DEVwarning(Sstr+"duplicate input of"+(Category+ReadInput::categoryTerminator+Name)
                            +"at lines"+(first+1)+"and"+(first+second+1)+"(first will be used)",1);
#ifndef _DEVELOP_
    ABORT(Sstr+"duplicate input of"+(Category+ReadInput::categoryTerminator+Name)
          +"at lines"+(first+1)+"and"+(first+second+1)+" --- remove/comment to continue");
#endif
    return true;
}

/// return input string
string ReadInput::readValue(const string Category, const string Name, const string Default, const string Docu, unsigned int Line, string Flag, string Allow) {
    // we need to up-cast, just in case...
    if(auto inp=dynamic_cast<ReadInputList*>(this))return inp->ReadInputList::readValue(Category,Name,Default,Docu,Line,Flag,{"allow all"});
    else                                           return      ReadInput::readItem(Category,Name,Default,Docu,Line,Flag,Allow)._value;
}

// locate the exact match of item name (between commas)
size_t ReadInput::posExactName(std::string Line, std::string Name){
    size_t pos(0);
    while((pos=Line.find(Name,pos))!=std::string::npos){
        size_t len=min(Line.find_first_of(itemSeparator+"=",pos),Line.length())-pos;
        if(tools::cropString(Line.substr(pos,len))==Name)break;
        pos+=len;
    }
    return pos;
}

InputItem ReadInput::dummyItem("DUMMY","DUMMY",0,"DUMMY","dummy item",ReadInput::noItem,"");
InputItem & ReadInput::readItem(const string Category, const string Name, const string Default, const string Docu, unsigned int Line, string Flag, string Allow) {

    if(dynamic_cast<ReadInputList*>(this)){
        for(auto l: allLines){
            if(l.find(Category+ReadInput::categoryTerminator+Name+"["+tools::str(Line)+"]=")==0){
                inputTable.push_back(InputItem(Category,Name,Line,Flag,Docu,Default,Allow));
                l=l.substr(l.find("]=")+2);
                l=l.substr(0,l.find('\\'));
                inputTable.back()._value=tools::cropString(l);
                return inputTable.back();
            }
        }
        error("input item not in "+_outDir+ReadInput::inputList,1);
    }
    duplicateInput(allLines,Category,Name);
    string itemDefinition=_root+Category+Name+tools::str(Line)+Docu;


    std::string inputString="UNDEFINED";
    string returnString=Default;
    inputTable.push_back(InputItem(Category, Name, Line, Flag, Docu, Default,Allow));

    if (Category[0] == ' ') error("Category must not start with a blank",inputTable.size()-1);
    if (string::npos != Category.find(categoryTerminator)) error("Do not include a colon in your category argument",inputTable.size()-1);
    if (string::npos != Name.find(categoryTerminator)) error("Do not include a colon in your input",inputTable.size()-1);
    if (string::npos != comments.find(Category[0])) error("Do not use comment characters '"+comments+"'' in your category argument",inputTable.size()-1);
    if (string::npos != comments.find(Name[0])) error("Do not use comment characters '"+comments+"'' in your input name",inputTable.size()-1);

    // postion of name in category line
    size_t namPos=std::string::npos;

    // command line arguments prevail, overwrite everything
    if(inCommandLine(Flag,Category+categoryTerminator+Name,Line,returnString)){
        inputString=returnString;
        goto Return;
    }


    // search for standard input item
    unsigned int n;
    for (n = 0; n < allLines.size(); n++){
        int column = -1;
        size_t pos = allLines.at(n).find(Category + categoryTerminator);
        if(comments.find(allLines.at(n)[0])!=string::npos)pos=string::npos; // commented out
        inputTable.back().lineFile = n;


        if(pos!=std::string::npos){

            if (pos != 0 and allLines.at(n).substr(pos-1)==" " and comments.find(allLines.at(n)[0])==string::npos)
                error("Found category name but not at first column; admissible comment characters: "+comments,inputTable.size()-1);

            // only detect presence of category
            if(Default==noItem and Name==anyName){
                returnString="FOUND:"+Category;
                inputString=returnString;
                goto Return;
            }

            // get for exact name in category
            namPos = posExactName(allLines.at(n),Name);

            // only category needed
            if(Name==anyName)namPos=0;

            if (namPos != string::npos){
                if(Default==flagOnly)error("specify input item exclusively through command line flag \"-"+Flag+"\"",inputTable.size()-1);

                // only detect presence of category and name
                if(Default==noItem){
                    returnString="FOUND:"+Category+Name;
                    inputString=returnString;
                    goto Return;
                }

                // value may be given as item=value
                size_t valPos=allLines.at(n).find("=", namPos);
                size_t itemEnd=findFirstOutsideBrackets(allLines.at(n),itemSeparator,leftBracket,rightBracket,namPos);
                if(valPos<itemEnd){
                    if(Line>1){
                        // treat name=value like a single input line, i.e. return default for any further line
                        if(allLines.size()>n+1 and allLines.at(n+1)!=""
                                and not isCategoryLine(allLines.at(n+1))
                                and comments.find(allLines.at(n+1)[0])==std::string::npos)
                            ABORT(Category+categoryTerminator+Name+" is multi-line input, cannot be specified as 'name=value'");
                        goto Return;
                    }
                    returnString=allLines.at(n).substr(valPos+1,itemEnd-valPos-1);
                    inputString=returnString;
                    goto Return;
                }
                if (n + Line >= allLines.size()){
                    inputString=ReadInput::notFound;
                    goto Return;
                }

                if(tools::subStringCount(allLines[n],",",leftBracket,rightBracket)<tools::subStringCount(allLines[n+Line],",",leftBracket,rightBracket))
                    ABORT("\nmore comma-separated inputs than items:\n   "+allLines[n]+"\n   "+allLines[n+Line]
                            +"\n Maybe you wanted to separate by blanks? Else enclose in single quotation marks:\n   \'"+
                            allLines[n+Line]+"\'");

                column = 0;
                while (string::npos != allLines.at(n).rfind(itemSeparator, namPos - 1)){
                    namPos = allLines.at(n).rfind(itemSeparator, namPos - 1);
                    column++;
                }

                // extract value string from line
                inputTable.back().lineFile = n + Line; // expected line for item
                returnString = allLines.at(n + Line);
                if(returnString.find_first_of(forbidden)!=string::npos)
                    error("Do not use any of the following characters in input value: "+forbidden+" (use single quote for quoting)",inputTable.size()-1);

                // locate columns
                for (int c = 0; c < column; c++){ // get the right column
                    size_t pos = findFirstOutsideBrackets(returnString,itemSeparator,leftBracket,rightBracket);
                    if (pos == string::npos){
                        returnString=Default;
                        inputString="EMPTY";
                        if (Default == noDefault) error("too few data items in line",inputTable.size()-1);
                    }
                    else {
                        returnString=returnString.substr(pos + 1);
                    }
                }
                // truncate other columns
                pos = findFirstOutsideBrackets(returnString,itemSeparator,leftBracket,rightBracket);
                if (pos!=string::npos)returnString=returnString.substr(0, pos);

                // remove trailing and leading blanks
                returnString=tools::cropString(returnString);
                if (returnString=="")returnString=Default;

            }
        }
        if (column > 0)break; // found the input
        inputString=notFound;
    }
    if(namPos==std::string::npos)inputString=notFound;

    if(inputTable.size()==0)ABORT("nothing read and no default given - probably using ReadInput::exclude without definition of the parameters");
    if (returnString==noDefault and not noInputFile())
        error("Input mandatory - category and name not found",inputTable.size()-1);

Return:
    removeQuotesAndBlanks(returnString);
    removeQuotesAndBlanks(inputString);
    if(Default==noItem){
        inputTable.pop_back();
        dummyItem._value=returnString;
        dummyItem._inpVal=inputString;
        return dummyItem;
    }
    else {
        if(readList.count(itemDefinition)!=0 and readList[itemDefinition]!=returnString)
            DEVABORT("input value of "+Category+":"+Name+"=\""+readList[itemDefinition]+"\" has changed to \""+returnString
                     +"\" -- must not change input values after start of program");
        readList[itemDefinition]=returnString;
        inputTable.back()._inpVal=returnString==Default?inputString:returnString;
        inputTable.back()._value=returnString;
        if(not inputTable.back().admissible())
            error("\n"+inputTable.back()._category+": "
                  +inputTable.back()._name+" ...\""+returnString+"\" not admissible"+
                  "\nAllowed inputs:\n"+strAdmissibleInputs(inputTable.back()._docu),inputTable.size()-1);

        if(inputTable.back().value()==flagOnly)inputTable.back()._value="false";
        if(inputTable.back().value()==flagFound)inputTable.back()._value="true";

        InputItem & it=findMatching(inputTable.back());
        if(&it!=&inputTable.back() or it.docu()==doNotAddToItemTable)inputTable.pop_back();
        return it;
    }
}

void ReadInput::exclude(string Cat1, string Cat2, string Nam1, string Nam2){
    if(dynamic_cast<ReadInputList*>(this))return; // machine generated file is always legal
    if(inputTable.size()==0)DEVABORT("must not use ReadInput::exclude before reading: "+Cat1+" <-> "+Cat2);
    // get flag names
    string flag1,flag2;
    for(InputItem i: inputTable){
        if(i._category==Cat1 and i._name==Nam1)flag1=i._flag;
        if(i._category==Cat2 and i._name==Nam2)flag2=i._flag;
    }

    string found1=readValue(Cat1,Nam1,noItem,"dummy",1,flag1,"");
    string found2=readValue(Cat2,Nam2,noItem,"dummy",1,flag2,"");
    if(     (found1!=noItem or found1.find("FOUND")!=string::npos) and
            (found2!=noItem or found2.find("FOUND")!=string::npos) ){
        string str1=Cat1+":"+Nam1; if(flag1!="")str1+=" [-"+flag1+"]";
        string str2=Cat2+":"+Nam2; if(flag2!="")str2+=" [-"+flag2+"]";
        if(Nam1!=anyName or Nam2!=anyName)
            ABORT("Input items mutually exclusive: "+str1+" - "+str2+", specify only one")
                    else
                    ABORT("Input categories mutually exclusive: "+Cat1+" - "+Cat2+", specify only one")
    }
}

bool ReadInput::hasCommandLineFlag(std::string Name){
    for(auto c: com)
        if(string(c).find("-"+Name)!=string::npos)return true;
    return false;
}

bool ReadInput::found(string Category, string Name, string Flag){
    string dum;
    if(inCommandLine(Flag,Category+categoryTerminator+Name,0,dum))return true;
    if(Flag==anyName and inCommandLine(Category+categoryTerminator+Name,Category+categoryTerminator+Name,0,dum))return true;
    if(readValue(Category,Name,noItem,"dummy",0,"","")=="FOUND:"+Category+Name)return true;
    return readValue(Category,Name,noItem,"dummy",0,"","")=="FOUND:"+Category;
}

bool ReadInput::isCategoryLine(string Line){
    if(Line.substr(0,Line.find_first_of(comments))=="")return false;
    if(tools::findFirstOutsideBrackets(Line,categoryTerminator,leftBracket,rightBracket)==string::npos)return false;

    // specifiction with "category[...]:" is allowed
    if(Line.find(categoryTerminator)>Line.find_first_of(notInCategory))return false;
    return true;
}

std::vector<std::string> ReadInput::categoryWithSpecification(std::string Category){
    std::vector<std::string> res;
    for(auto l: allLines){
        if(l.find(Category+"[")==0 and l[l.find("]")+1]==':')res.push_back(l.substr(0,l.find(":")));
    }
    return res;
}

std::set<std::string> ReadInput::namesWithSpecification(std::string Category, std::string Name){
    std::set<std::string> res;
    for(auto l: allLines){
        if(isCategoryLine(l) and l.find(Category+categoryTerminator)==0){
            std::string nameSpec;
            for(size_t pos=0;(pos=l.find(Name+"[",pos))!=std::string::npos;pos+=nameSpec.size()){
                auto beg=l.find("[",pos)+1,end=l.find("]",pos);
                std::string spec=tools::cropString({l.begin()+beg,l.begin()+end});
                nameSpec=Name+"["+spec+"]";
                if(not found(Category,nameSpec)){
                    ABORT("incorrect syntax of name-specification of "+Category+categoryTerminator+nameSpec
                          +"\nin Line\n"+l+"\ndo not use blanks inside [...]");
                }
                res.insert(spec);
            }
        }
    }
    return res;
}

bool ReadInput::endCategory(string Category, int Line) {

    if(not found(Category))return true;

    for(size_t n=0;n<allLines.size();n++){
        if(allLines[n].find(Category+categoryTerminator)!=string::npos
                and isCategoryLine(allLines[n]))
        {
            if(size_t(n+Line)>=allLines.size())return true;
            // category terminates by blank or next category
            for(int ll=1;ll<=Line;ll++){
                if(tools::cropString(allLines[n+ll])=="" or isCategoryLine(allLines[n+ll]))return true;
            }
            return false;
        }
    }
    return true;
}

bool ReadInput::inCommandLine(string Flag, string DefaultFlag, unsigned int Line, string & valString){

    // no flag defined for item
    if(Flag=="")return false;

    auto itFlag=com.begin();
    while(itFlag!=com.end() and ((*itFlag+"=").find("-"+Flag+"=")!=0))itFlag++;
    if(itFlag==com.end()){
        for(itFlag=com.begin();itFlag!=com.end();itFlag++)
            if(itFlag->find("-"+DefaultFlag+"=")==0)
                ABORT("use abbreviated flag -"+Flag+" for "+*itFlag);
        return false;
    }

    int flagPos=itFlag-com.begin();

    // locate value(s) after flag
    size_t valPos=com[flagPos].find("=");
    if(valPos==string::npos)valString=flagFound;
    else valString=com[flagPos].substr(valPos+1);
    if(valString=="")valString=flagFound;

    // advance to Line'th value in input list
    if(Line>1)DEVABORT("multi-line input from command line not functional");
    if(Line>1 and tools::findFirstOutsideBrackets(valString,":","[(","])")==string::npos)
        ABORT(Str("multi-line item at line ","")+Line+
              ",\n command line input as -"+Flag+"="+string(':',Line-1)+"value:.. (after "+(Line-1)
              +"th \':\')\n input string is: "+valString);
    valPos=0;
    for(unsigned int l=1;l<Line;l++){
        valPos=tools::findFirstOutsideBrackets(valString,":","[(","])",valPos+1);
        // Line'th item not specified, use standard input
        if(valPos==string::npos)return false;
    }

    // cut out input value, remove leading and trailing blanks
    valString=tools::cropString(valString.substr(valPos,tools::findFirstOutsideBrackets(valString,":","[(","])",valPos+1)));
    return true;
}

void ReadInput::inputMacros(){
    inpLines=allLines; // keep a copy of the original input lines

    // check macro definitions for consistency
    std::map<string,int> allMacros;
    for(size_t k=0;k<allLines.size();k++){
        size_t def=allLines[k].find(macro);
        if(def==string::npos)continue; // not a macro line
        string macr=tools::cropString(allLines[k].substr(macro.length()));
        if(allMacros.find(macr)!=allMacros.end())
            ABORT(Str("multiple definition: #define")+macr+"in line"+(k+1)+"first found in line"+(allMacros[macr]+1));
        if(def!=0)ABORT("found macro definition, but not starting in first column\n"+allLines[k]);
        allMacros[macr]=k;
        if(macr.find(" ")==string::npos)
            ABORT("ill-formed macro definition at line "+tools::str(k)
                  +"\n"+tools::cropString(allLines[k])
                  +"\nUse blank to separate macro name from value");
    }

    for(unsigned int k=0;k<allLines.size();k++){
        size_t def=allLines[k].find(macro);
        if(def!=string::npos){
            string macr=tools::cropString(allLines[k].substr(macro.length()));
            string subs=tools::cropString(macr.substr(macr.find(" ")));
            macr=macr.substr(0,macr.find(" "));
            for(unsigned int l=k+1;l<allLines.size();l++){
                // search allLines for macro and recursively substitute if found found
                size_t pos;
                while(string::npos!=(pos=allLines[l].find(macr)))
                    allLines[l]=allLines[l].substr(0,pos)+subs+allLines[l].substr(pos+macr.length());
            }
        }
    }
}


string ReadInput::lineAt(string Category, int Line){
    if(endCategory(Category,Line))return notFound;
    int cnt=Line;
    for(std::string l: allLines){
        if(cnt<Line or (cnt==Line and isCategoryLine(l) and l.find(Category+categoryTerminator)!=string::npos))cnt--;
        if(cnt==-1)return l;
    }
    return notFound;
}
std::string ReadInput::outputFile(){return output()+PrintOutput::outExtension;}

string ReadInput::outputTopDir() const {
    if(_outDir!="NONE")return _outDir.back()=='/'?_outDir:_outDir+"/";

    if(_root=="")return "";

    if(MPIwrapper::Rank()==int(mpiMaster)){
        // try create directory
#ifdef _WIN32
        if(not folder::exists(root))
            if (not folder::create(root)) error("could not create output main directory "+(root));
#else
        if (not folder::exists(_root))
            if (not folder::create(_root)) error("could not create output main directory "+(_root)+"\n may exist but does not contain 0000 subfolder",inputTable.size());
#endif
        _outDir=nextRunDir(_root);

        // create directory
        if (not folder::exists(_outDir))
            if (not folder::create(_outDir))error("could not create output main directory " + (_outDir),inputTable.size());

        // append slash to the root file
        if(_outDir.rfind("/")!=_outDir.length()-1)_outDir+="/";

        // write copy of input into output directory
        ofstream inpc((_outDir+inputCopy).c_str());
        for (unsigned int n=0;n<inpLines.size();n++)
        {
            inpc<<inpLines.at(n)<<endl;
        }
        inpc.close();
    }
    MPIwrapper::Bcast(_outDir,mpiMaster);

    //    if(not MPIwrapper::isMaster())outDir+=tools::str(MPIwrapper::Rank())+"_";
    return _outDir;
}

void ReadInput::outSubdir(std::string Name){
    // create path+Name of subdirectory
    subDir = Name+"/";
    std::string _outSubdir = outputTopDir()+Name;

    // outDir = "tutorial/20Helium2d/0123/"
    // root   = "tutorial/20Helium2d"

    if(MPIwrapper::Rank()==int(mpiMaster)){
        // try create directory (check whether root exists)
#ifdef _WIN32
        if(not folder::exists(outputTopDir()))
            if (not folder::create(outputTopDir())) error("could not create output subdirectory "+(outputTopDir()));
#else
        if (not folder::exists(outputTopDir()))
            if (not folder::create(outputTopDir())) error("could not create output subdirectory "+(outputTopDir()),inputTable.size());
#endif

        // create directory
        if (not folder::exists(_outSubdir))
            if (not folder::create(_outSubdir))error("could not create output main directory " + (_outSubdir),inputTable.size());

        // append slash to the root file
        if(_outSubdir.rfind("/")!=_outSubdir.length()-1)_outSubdir+="/";

    }
    MPIwrapper::Bcast(_outSubdir,mpiMaster);

    if(not MPIwrapper::isMaster())_outSubdir+=tools::str(MPIwrapper::Rank())+"_";
}

void ReadInput::show() const {
    cout << "\n---- beginning of file " << (_root) << ReadInput::inputExtension << " --------------------------------------------------------" << endl;
    for (unsigned int n = 0; n<allLines.size(); n++)cout << allLines.at(n) << endl;
    cout << "---- end of file --------------------------------------------------------" << endl;
}

void ReadInput::showHelp(std::vector<string> CommandLine, std::string Category) {
    // this really could go into its own little "help()"
    if(MPIwrapper::isMaster()){
        std::string category(Category);
        if(category=="" and CommandLine.size()>1){
            for (std::string c: CommandLine){
                if(c.find("-h")==0 or c.find("--help")==0){
                    if(c.find("=")!=string::npos){
                        category=c.substr(c.find("=")+1);
                    }
                }
            }
        }

        // print help message and terminate
        string line;
        if(CommandLine.size()){
            ifstream stream((CommandLine[0]+docExtension).c_str(),std::ios::in);
            if(stream.is_open()){
                bool show=false;
                if(category!=""){
                    cout<<endl;
                    while(getline(stream,line)){
                        // start printing at summary
                        show=show or line.find(category+ReadInput::categoryTerminator)==0;
                        if(show){
                            if(line=="" or line.find(markerSummary)!=std::string::npos)break;
                            std::cout<<line<<std::endl;
                        }
                    }
                }
                else {
                    while(getline(stream,line)){
                        // start printing at summary
                        show=show or line.find(markerSummary)!=std::string::npos;
                        if(show)std::cout<<line<<std::endl;
                    }
                    cout<<"\n(--- end "+CommandLine[0]+docExtension+" ---)"<<endl;
                }
            }
        }
        std::cout<<std::endl;
        PrintOutput::warning(
                    "No input file specified, \nusage: > "+
                    CommandLine[0]+" InputFile [-flag1=val1 -flag etc.]\n       > "+
                CommandLine[0]+" -h=Categ ... specific help for Categ\n");
    }
    MPIwrapper::Finalize();
    if(Category=="")exit(1);

}

// structured error message
void ReadInput::error(string Message, unsigned int Item,unsigned int Item2) const {
    cout << "\n+++ Input error: " + Message + "\nFile: " << _root << ReadInput::inputExtension << endl;
    if (inputTable.size()>Item){
        cout << "Category and Name: " + inputTable[Item]._category + ": " + inputTable[Item]._name
                +" ["+inputTable[Item]._docu+"] -"+inputTable[Item]._flag<< endl;
        if(Item2<inputTable.size())cout << "Category and Name: " + inputTable[Item2]._category + ": " + inputTable[Item2]._name
                                           +" ["+inputTable[Item2]._docu+"] -"+inputTable[Item2]._flag<< endl;
        if(allLines.size()>inputTable[Item].lineFile)
            cout << "Line " << inputTable[Item].lineFile << ": " << allLines.at(inputTable[Item].lineFile)<<" "<< endl;
    }
    if(Item<inputTable.size())showHelp(commandLine,inputTable[Item].category());
    ABORT(Message);
}

// conversion to bool: must be either true or false
bool ReadInput::string_to_bool(const string &Text){
    if (Text != "true"){
        if (Text != "false")error("bool string must be true or false, is: \"" + Text+"\"",inputTable.size()-1);
        return false;
    }
    return true;
}

// conversion to double and int
int    ReadInput::string_to_int(const string &Text){ return atoi(Text.c_str()); }
double ReadInput::string_to_double(const string &Text){
    if(Text.find("inf")!=string::npos){
        if(Text[0]=='-')return -DBL_MAX;
        else            return DBL_MAX;
    }
    char* pEnd; return strtod(Text.c_str(),&pEnd); }

void ReadInput::setUnits(string Units){unitSystem=Units;}

void ReadInput::parseInput(std::string Line, std::string &Category, std::vector<std::string> &Name){
    // split input line into categories and names
    Category="";
    Name.clear();
    if (Line[0]==' ')return; // not a category line

    // ignore categoryTerminator after other non-alphanumeric characters
    if(Line.find(categoryTerminator)>=Line.find_first_of("!([<|{,;/#"))return;

    // ignore catagoryTerminator in between quotes or brackets
    string::size_type sep = findFirstOutsideBrackets(Line,categoryTerminator,quote+"([<|{",quote+")]>>}"),assign;

    if (string::npos == sep)return; // not a category line
    Category = Line.substr(0, sep);
    Line = tools::lcropString(Line.substr(sep + 1));
    while (Line != ""){
        Name.push_back(Line);
        sep=findFirstOutsideBrackets(Line,itemSeparator,leftBracket,rightBracket);
        assign=Line.find("=");
        string::size_type cut=min(sep,assign);
        Name.back() = tools::rcropString(Name.back().substr(0, cut));
        if (sep == string::npos)break; // no further names
        Line = tools::lcropString(Line.substr(sep + 1));
    }
}

static bool inputItemLess(const InputItem & a, const InputItem & b){
    if (a.category() != b.category())return a.category() < b.category();
    return a.name() < b.name();
}

static std::map<std::string,std::string> _texdocu;
static std::map<std::string,std::string> _categorySort;

void ReadInput::texdocuCategoryAdd(std::string Category, std::string Sort, std::string TexString,std::string Tutorials){
    LatexDocu::categoryAdd(Category,Sort,TexString,Tutorials);
}

void ReadInput::writeLinp() const {
    if(not _writeLinp)return;
    string cat="",name;
    std::ofstream list((outputTopDir()+inputList).c_str());

    list<<"--- DO NOT EDIT: InputItem's as actually appearing in code -----"<<endl;

    for (unsigned int n = 0; n < inputTable.size(); n++){
        if (cat != inputTable.at(n)._category){
            // new category - reset name
            cat = inputTable.at(n)._category;
            name="";
        }
        std::string mutab=inputTable.at(n)._docu.find("(mutable)")!=std::string::npos?UiInputItem::markUp("mutable",""):"";
        list<<inputTable.at(n).strMarkup()<<mutab<<endl;
        if(inputTable[n]._name==name)continue; // document each name only once
        name=inputTable[n]._name;
    }
    list<<std::flush;
}

void ReadInput::writeDoc(ostream *Doc, std::vector<string> AddLines){

    LatexDocu::write(defaultLaTeXdir,inputTable);
    if(not MPIwrapper::isMaster(MPIwrapper::worldCommunicator()))return;

    ostream *doc=Doc;
    std::vector<std::string> previousLinp;
    if(folder::exists(outputTopDir()+inputList)){
        std::string line;
        ifstream prevLinp(outputTopDir()+inputList);
        while(getline(prevLinp,line))previousLinp.push_back(line);
    }

    if(doc==0)doc=new ofstream((mainProgram+docExtension).c_str());
    *doc<<"\n ---------------------\n ! Admissible inputs ! \n ---------------------\n"<<endl;
    *doc<<"Category:\nName (-CommandLineFlag) [defaultValue]\n"<<endl;
    if(com.size()==0)*doc<<"\n!!!! command line flags DISABLED - supply argc,argv from main(int argc,char* argv[]) when creating input !!!\n";
    string cat="",name;
    *doc<<"you can define macros for the input file in the form\n"+macro+" macroString replacementString"<<endl;
    *doc<<"\nA list of all inputs will be written to output_directory/"<<inputList<<endl;
    if(_root==flagOnly)*doc<<"\n  !!! --- ONLY FLAG INPUT ALLOWED --- !!!"<<endl;
    stable_sort(inputTable.begin(), inputTable.end(), inputItemLess);

    for (unsigned int n = 0; n < inputTable.size(); n++){
        if (cat != inputTable.at(n)._category){
            cat = inputTable.at(n)._category;
            *doc << endl << cat << ":" << endl;
            name="";
        }
        if(inputTable[n]._name==name)continue; // document each name only once
        name=inputTable[n]._name;
        *doc<<inputTable[n].docuLine()<<std::endl;
    }

    // write list of all inputs for machine read
    writeLinp();


    if(AddLines.size()>0){
        *doc<<"\n------- Additional info ----------------------------";
        for(auto l: AddLines)*doc<<"\n "+l;
        *doc<<std::endl;
    }

    *doc<<"\n------------------ Summary of input categories ------------------------";
    cat="";
    for (unsigned int n = 0; n < inputTable.size(); n++){
        if(cat!=inputTable.at(n)._category){
            cat=inputTable.at(n)._category;
            *doc<<"\n"<<cat<<categoryTerminator<<" ";
        } else {
            if(inputTable.at(n)._name!=name)*doc<<", ";
        }
        if(inputTable.at(n)._name==name)continue;
        name=inputTable.at(n)._name;
        *doc<<name;
    }
    *doc<<endl;
    PrintOutput::message("Description of input on file "+mainProgram+docExtension,&cout,true,1);
    if(Doc==0)delete doc;

    // check change of input
    if(previousLinp.size()>0){
        std::string line;
        ifstream currLinp(outputTopDir()+inputList);
        int k=0;
        std::vector<std::string> newLines;
        while(getline(currLinp,line)){
            if(std::find(previousLinp.begin(),previousLinp.end(),line)==previousLinp.end())newLines.push_back(line);
#ifndef _DEVELOP_
            if(line.find("DEBUG")!=std::string::npos and newLines.size()>0)newLines.pop_back(); // do not confuse standard user with DEBUG inputs
#endif
            k++;
        }
        if(k<int(previousLinp.size()) or newLines.size()){
            bool mute=true;
            for(auto l: newLines)mute=mute and
                    (  l.find("(mutable)")!=std::string::npos
                    or l.find("\\mutable{")!=std::string::npos // new markup style
                    or l.find("DEBUG")!=std::string::npos);
            if(not mute and _writeLinp){
                std::string prevLinp=tools::newFile(outputTopDir()+inputList);
                std::ofstream prevStream(prevLinp);
                for(auto l: previousLinp)prevStream<<l<<"\n";

                PrintOutput::title("inputs changed");
                PrintOutput::warning("inputs not defined as mutable - previous saved in "+prevLinp);
                if(newLines.size())PrintOutput::DEVmessage("suppress linp-file copy by adding (mutabel) to read docu string");
            }
            else
                PrintOutput::title("mutable inputs changed");

            PrintOutput::paragraph();
            for(auto l: newLines)PrintOutput::subTitle(" - "+l.substr(0,l.find("\\")));

        }
        else {
            if(k<int(previousLinp.size()))
                PrintOutput::message(Sstr+"fewer inputs: this="+k+"previous="+previousLinp.size());
        }
    }
}

string ReadInput::allowedString(string Value, InputItem & It){
    if(It._allow.length()==0)return Value;
    vector<string> allVal=tools::splitString(It._allow,',');

    allVal.push_back(It._defVal);
    for(unsigned int k=0;k<allVal.size();k++)
        if(tools::cropString(allVal[k])==Value){
            return Value;
        }
    ABORT("illegal input \""+Value+"\" allowed values: "+It._defVal+"(Default),"+It._allow);
}

static double value_from_alg(std::string Value){
    Algebra alg(Value);
    if(not alg.isAlgebra())ABORT("not a number or algebraic expression: "+Value);
    return alg.val(0.).real();
}

double ReadInput::allowedNumber(string Value, InputItem& It){

    double value=value_from_alg(Value);

    if(It._allow.length()==0)return value;
    if(It.admissible())return value;

    vector<string> allVal=tools::splitString(inputTable.back()._allow,',',"[","]");
    allVal.push_back(inputTable.back()._defVal);

    for(unsigned int k=0;k<allVal.size();k++){
        string crop=tools::cropString(allVal[k]);
        if(crop.find("[")==0){
            vector<string> ab=tools::splitString(tools::stringInBetween(crop,"[","]"),',');
            if(ab.size()!=2)ABORT("specify range of allowed values as \"[a,b]\", found: "+crop);
            if(tools::string_to_double(ab[0])<=value and value<=tools::string_to_double(ab[1]))return value;
        } else {
            if(crop==Value)return value;
        }
    }
    error("illegal input \""+tools::str(Value)+"\" allowed values: "+It._defVal+" (Default), "+It._allow,0);
    return 0.;
}

void ReadInput::finish(){
    // check sanity of input Categories and Names

    if(not MPIwrapper::isMaster())return;
    writeLinp();// write all inputs for easy parsing (hard human reading)
    writeDoc(); // generate current input documentation

    // if no command line is supplied, flags cannot be used
    if(com.size()==0)
        for (unsigned int k=0;k<inputTable.size();k++)
            if(inputTable[k].flag()!=inputTable[k]._category+categoryTerminator+inputTable[k]._name){
                cout<<"\n ERROR: flag input enabled, but command line not supplied \n     supply narg,argv[] when creating input"<<endl;
                abort();
            }

    // input items must be unique
    for(int n=inputTable.size()-1;n>=0;n--){
        for(int m=0;m<n;m++){
            if(inputTable[m]._category==inputTable[n]._category and inputTable[m]._name==inputTable[n]._name){
                // remaining parameters must match
                if(inputTable[m]._flag!=inputTable[n]._flag
                        or inputTable[m]._docu!=inputTable[n]._docu)
                    error("duplicate use of input Category and Name",n,m);
            }

            // flags must be unique through all inputs
            if(inputTable[m]._flag==inputTable[n]._flag and inputTable[m]._flag!=""){
                if(inputTable[m]._category!=inputTable[n]._category
                        or inputTable[m]._name!=inputTable[n]._name)
                    error("duplicate use of command line flag '-"+inputTable[m]._flag+"' ",n);
            }

        }
    }
    // check command line for legal flags
    bool found;
    for (int flagPos=0;flagPos<(int)com.size();flagPos++){
        if(com[flagPos][0]!='-')continue;
        string flag=tools::cropString(com[flagPos].substr(1,com[flagPos].find("=")-1));
        found=false;
        for (unsigned int k=0;k<inputTable.size();k++){
            found=inputTable[k]._flag==flag;
            if(found)break;
        }

        if (not found and flag.find("DEBUG")!=0  and flag.find("AMP")!=0 and flag!="h"){
            if(flag=="help" or flag=="-help" or flag=="h"){
                writeDoc(&cout);
                exit(0);
            }

            else {
                cout << "undefined flag '"+flag+"' on command line: "+com[flagPos]<< endl;
                cout << "check file " <<mainProgram<<docExtension << " for admissible inputs" << endl;
                PrintOutput::message("check file "+mainProgram+docExtension+" for admissible inputs");
                PrintOutput::DEVmessage("ReadInput::finish() may have been called before input was actually finished");
                abort();
            }
        }
    }

    string cat;
    vector<string> name;
    for (unsigned int n = 0; n < allLines.size(); n++){
        // category string from line
        parseInput(allLines.at(n), cat, name);
        found = false;
        for (unsigned int l = 0; l < name.size(); l++){
            for (const InputItem & kk: inputTable){
                found = kk._category == cat and kk._name == name[l];
                if(found)break;
            }

            if (not found and comments.find(allLines.at(n)[0])==string::npos){
                //                showHelp();
                cout << "\n"+allLines.at(n) << endl;
                cout << "\ncategory or name '"<<cat<<ReadInput::categoryTerminator<<" "<< name[l]<<endl;
                PrintOutput::message("check file "+mainProgram+docExtension+" for admissible inputs");
                PrintOutput::DEVmessage("ReadInput::finish() may have been called before input was actually finished");
                ABORT("suspicious input line "+tools::str(n+1));
            }
        }
    }
    if(noInputFile()){
        PrintOutput::message("called without input file - see "+mainProgram+docExtension+" for possible inputs or run \""+mainProgram+" -help\"");
        exit(0);
    }
}

void ReadInput::obsolete(string Category, string Name, string Message, unsigned int Line, string Flag)
{
    if(readValue(Category, Name, "OBSOLETE", "...OBSOLETE...", Line, Flag,"")!="OBSOLETE"){
        //        error("OBSOLETE INPUT",inputTable.size()-1);
        PrintOutput::title("OBSOLETE INPUT");
        PrintOutput::message(Category+categoryTerminator+" "+Name+(Flag==""?"":" (flag "+Flag+")"));
        PrintOutput::message(Message);
        PrintOutput::paragraph();
        exit(0);
    }
}

void ReadInput::flagError(string Flag) const
{
    error("specify value -"+Flag+"=???, without blanks around \"=\" ",inputTable.size()-1);
}

// overload conversion to final values
InputItem & ReadInput::read(string Category, string Name, string & Value,
                            string Default, string Docu, unsigned int Line, string Flag, string Allow)
{
    InputItem & it=readItem(Category, Name, Default, Docu, Line, Flag, Allow);
    string str(it._value);
    if(str==flagFound)flagError(Flag);
    Value=str;
    Value=tools::cropString(Value);
    Value=tools::unquote(Value);
    allowedString(ReadInputRange::low(Value),it);
    return it;
}
InputItem & ReadInput::read(std::string Category, std::string Name, int & Value,
                            std::string Default, std::string Docu, unsigned int Line, string Flag, string Allow)
{
    InputItem & it=readItem(Category, Name, Default, Docu, Line, Flag, Allow);
    string str(it._value);
    if (str==flagFound)flagError(Flag);
    Value = (int) allowedNumber(ReadInputRange::low(str),it);
    return it;
}
InputItem & ReadInput::read(std::string Category, std::string Name, unsigned int & Value,
                            std::string Default, std::string Docu, unsigned int Line, string Flag, string Allow)
{
    InputItem & it=readItem(Category, Name, Default, Docu, Line, Flag, Allow);
    string str(it._value);
    if(str==flagFound)flagError(Flag);
    Value = (int) allowedNumber(ReadInputRange::low(str),it);
    return it;
}

InputItem & ReadInput::read(std::string Category, std::string Name, double & Value,
                            std::string Default, std::string Docu, unsigned int Line, string Flag,string Allow)
{
    InputItem & it=readItem(Category, Name, Default, Docu, Line, Flag, Allow);
    string str(it._value);
    str=ReadInputRange::low(str);
    if(str==flagFound)flagError(Flag);
    Value = string_to_double(ReadInputRange::low(str));

    // check wether unitSystem or dimensions are specified
    str=tools::cropString(str);
    size_t lastDig=str.find_first_not_of("+-.0123456789e*/()");

    // no unitSystem specified, return value "as is"
    if(isInfinity(str))return it;

    // check for algebra
    string alg=str.substr(0,str.rfind(' '));
    alg=alg.substr(0,alg.find('~')); // units may be attached by ~

    Algebra::failures="";
    if(not Algebra::isAlgebra(alg))
        ABORT("\""+Category+": "+Name+"\" input \""+str+
              "\"\nnot an algebraic expression:\n   "+Algebra::failures
              +"\nAvailable constants: "+Algebra::listConstants()
              +"\nNote: write algebraic expression WITHOUT blanks, if units are specified, separate by blank or ~ (useful on command line)"
              +"\n      e.g. 2*pi/3 m or 2*pi/3~m for \"meters\""
              );
    Value=real(Algebra(alg).val(1.));

    if(str.length()==alg.length()){
        allowedNumber(str,it);
        return it; // only algebra, not units
    }

    // make sure standard units are set
    Units::standardUnits();
    lastDig--;

    string inp=tools::cropString(str.substr(lastDig+1));
    if(inp[0]=='~')inp=inp.substr(1); // alternative to blank, attach units by ~
    string mess;
    if(toLower(inp).find(toLower("Infty").substr(0,3))!=string::npos)
        mess="\n\n ??? maybe you meant \"infty\", found: \""+inp+"\"\n";

    if(not tools::hasKey(Units::aka,inp))error("unknown unit \""+inp+"\"\navailable: "+tools::listMap(Units::aka)+mess,ReadInput::inputTable.size()-1);

    // check Name for output unitSystem and dimension
    std::string out="";
    if(Name.find("(")!=string::npos)out=tools::stringInBetween(Name,"(",")");
    if(tools::hasKey(Units::aka,out))out=Units::aka[out]; // known unit name, like nm, W/cm2 etc.

    if(out=="")out=unitSystem;
    else if(out.find(Units::sep)==string::npos)out=unitSystem+Units::sep+out;
    if(out.find("UNDEFINED")!=string::npos)
        ABORT("for using input unitSystem, set default unitSystem for setUnits(name)");
    Value=Units::convert(Value,inp,out);
    allowedNumber(tools::str(Value),it);

    return it;
}


InputItem & ReadInput::read(std::string Category, std::string Name, bool & Value,
                            std::string Default, std::string Docu, unsigned int Line, string Flag)
{
    InputItem & it=readItem(Category, Name, Default, Docu, Line, Flag,"");
    // trim possible markups
    it._value=tools::cropString(it.value().substr(0,it.value().find_first_of(" \\")));
    if(it._value==flagOnly) it._value="false";
    if(it._value==flagFound)it._value="true";

    Value = string_to_bool(it.value());
    return it;
}

// vectors of input
InputItem & ReadInput::read(std::string Category, std::string Name, vector<string> & Value,
                            std::string Default, std::string Docu, unsigned int Line, string Flag, string Allow){
    string line;
    InputItem & it=read(Category,Name,line,Default,Docu,Line,Flag,Allow);
    line=tools::cropString(line)+" ";
    while(line!=" "){
        Value.push_back(allowedString(line.substr(0,min(line.find(" "),line.length())),it));
        line=tools::cropString(line.substr(line.find(" ")))+" ";
    }
    return it;
}

InputItem & ReadInput::read(std::string Category, std::string Name, vector<double> & Value,
                            std::string Default, std::string Docu, unsigned int Line, string Flag, string Allow)
{
    vector<string> line;
    InputItem & it=read(Category,Name,line,Default,Docu,Line,Flag,Allow);
    for(unsigned int i=0;i<line.size();i++)Value.push_back(allowedNumber(line[i],it));
    return it;
}

InputItem & ReadInput::read(std::string Category, std::string Name, vector<int> & Value,
                            std::string Default, std::string Docu, unsigned int Line, string Flag, string Allow)
{
    vector<string> line;
    InputItem & it=read(Category,Name,line,Default,Docu,Line);
    for(unsigned int i=0;i<line.size();i++)Value.push_back((int) allowedNumber(line[i],it));
    return it;
}


/// consult this subroutine for models how to use the class
void ReadInput::Test(int NCom, char * Com[]){

    // create a dummy input file first
    string inputfile = "example.inp";
    std::ofstream o(inputfile.c_str());
    o << "Example input file - almost any format of comments can be written into it" << endl;
    o << "" << endl;
    o << "Below you find a categories with a few items" << endl;
    o << " - A Category starts with a non-blank item in column 1 and ends with a categoryTerminator (" << categoryTerminator << ")" << endl;
    o << " - Item names follow after the colon and are separated by commas" << endl;
    o << " - Blanks at the beginning and end of item values will be removed" << endl;
    o << " - Items can be retrieved in any order" << endl;
    o << " - Items not found will return the default value" << endl;
    o << "        if default value 'ReadInput::noDefault' is specified, an exception occurs" << endl;
    o << endl;
    o << "My other category: other item" << endl;
    o << " 0.2 " << endl << endl;
    o << endl;
    o << "the first item below can be overwritten by a flag -first=value"<<endl;
    o << "My input category: first item, third item, second item" << endl;
    o << "0.3  ,third item is a string,  24 " << endl;
    o << endl;
    o << "# you may also specify values in line (only for single line input)"<<endl;
    o << "In line value: first, second, third=17" << endl;
    o << ",,this will be ignored"<<endl;
    o << ",24"<<endl;
    o << endl;
    o << "# any input can be overwritten by command line flag"<<endl;
    o << "Overwrite: first, second" << endl;
    o << ",42."<<endl;
    o.close();

    // open input file and read all lines into buffer
    ReadInput in(inputfile,NCom,Com);

    // open input file for "main" input
    ReadInput::openMain(inputfile,NCom,Com);

    // read input items in any desired sequence
    int alt1, second;
    double alt2, first;
    string third;
    in.read("My other category", "other item", alt2, "137", "defaults must be specified as strings for all input types");
    in.read("My input category", "second item", second, "42", "not present in input, default suplemented");
    in.read("My input category", "first item", first, "27.", "this is a silly first item",1,"first");
    // admit overwriting by flag (last argument)
    in.read("My input category", "third item", third, ReadInput::noDefault, "third item must be present");

    in.read("My other category", "first item", alt1, "31415", "this is an alternative other item");

    // the following line leads to error: duplicates the flag "first"
    //in.read("My other category", "wrong item", alt1, "31415", "this is an alternative other item",1,"first");

    double inLine1,inLine2,inLine3,inLine22,over1,over2;
    in.read("In line value", "first",  inLine1, "101.", "this is an alternative other item");
    in.read("In line value", "second", inLine2, "101.", "this is an alternative other item");

    // the following line leads to error: duplicate input with different doc or default value
    //    in.read("In line value", "second",  inLine22, "101.", "this is an alt***ative other item");

    in.read("In line value", "third",  inLine3, "101.", "this is an alternative other item");
    in.read("In line value", "second",  inLine22, "101.", "this is an alternative other item",2);

    in.read("Overwrite", "first", over1, "101.", "dummy example for overwriting");
    in.read("Overwrite", "second", over2, "101.", "another dummy");


    cout << "========================================================" << endl << endl;
    in.show(); // show the content of in
    cout << "\nReadInput read return values:" << endl;
    cout << "\nMy input category:\n";
    cout << "   first: " << first << " (can be overwritten by -first=value)" << endl;
    cout << "  second: " << second << endl;
    cout << "   third: " << third << endl;
    cout << "\nMy other category:\n";
    cout << "    alt1: " << alt1 << " (default value)" << endl;
    cout << "    alt2: " << alt2 << endl;
    cout << "\nIn line value:\n";
    cout << " inLine1: " << inLine1 <<" (default value)"<<endl;
    cout << " inLine2: " << inLine2 <<endl;
    cout << " inLine3: " << inLine3 <<" (inline value, line 1 value ignored)"<< endl;
    cout << "inLine22: " << inLine22 <<" (second line in column 2)"<< endl;
    cout << "\nOverwrite:\n";
    cout << "    second: " << over2 <<" (try overwriting by -Overwrite:second=value)"<< endl;
    cout << "========================================================" << endl << endl;

    // check input for consistency (misprints etc)
    in.finish();
}

bool ReadInput::operator==(const ReadInput & Other) const {
    if(this==&Other)return true;
    if(allLines.size()!=Other.allLines.size())return false;
    for(int k=0;k<allLines.size();k++)
        if(allLines[k]!=Other.allLines[k])return false;
    return true;
}
