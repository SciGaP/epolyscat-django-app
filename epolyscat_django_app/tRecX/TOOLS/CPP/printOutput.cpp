// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "printOutput.h"
#include <sstream>
#include "tools.h"
#include "folder.h"
#include <float.h>

#include <fstream>
#include <iostream>
#include "multiStream.h"

#include <memory>
#include <outputTable.h>

#include "asciiFile.h"
#include "readInput.h"
#include "mpiWrapper.h"

using namespace std;

const string PrintOutput::outExtension="outf";
string PrintOutput::fileOut;
string PrintOutput::monMessage;
std::ofstream * PrintOutput::sOut=0;
std::ofstream * PrintOutput::sMon=0;
MultiStream * PrintOutput::out=new MultiStream();
unsigned int PrintOutput::lengthLine=80;
unsigned int PrintOutput::namesTab=20;
string PrintOutput::_indent="   ";
string PrintOutput::titleSeparator="=";
string PrintOutput::itemSeparator=", ";
string PrintOutput::nameSeparator=",";
string PrintOutput::_columnSeparator=" ";
string PrintOutput::namesEnd=":";

map<string,unsigned int> PrintOutput::warningCount;
map<string,unsigned int> PrintOutput::messageCount;

vector<string> PrintOutput::finalMessage;
string PrintOutput::level="full";

string PrintOutput::namesLine="",PrintOutput::valuesLine="";

// file-local variables
static string progressHeader="";
static string progressStatus="";
static vector<string> previousLevels;

bool PrintOutput::noScreenFlag() {
    bool noScreen;
    ReadInput::main.read("Print","noScreen",noScreen,ReadInput::flagOnly,"suppress output to screen",0,"noScreen");
    return noScreen;
}

void PrintOutput::set(string File){

    if(not MPIwrapper::isMaster())return;
    if(File.length()==0 or File=="")return;
    paragraph();
    paragraph();
    if(File=="cout" or File==outExtension){
        out->remove("Screen");
        out->add("Screen",&std::cout);
        fileOut="NONE";
    } else {
        if(folder::exists(File)){
            int i=0;
            std::string fileC(File);
            while(folder::exists(fileC))fileC=File+tools::str(i++);
            AsciiFile(File).copy(fileC,true);
            PrintOutput::message("existing PrintOutput file copied to "+fileC);
        }
        message("further output on "+File);
        fileOut=File;
        sOut=new ofstream();
        sOut->open(File.c_str());
        out->replace("File",sOut);
    }
}

void PrintOutput::outputLevel(string Level){
    if(Level!="full" and Level!="low" and Level!="off" and Level!="restore")
        ABORT("output level must be one of full,low,off,restore is: "+Level);
    if(Level=="restore"){
        if(previousLevels.size()>0){
            level=previousLevels.back(); // simple solution, later push/pop back
            previousLevels.pop_back();
        }
    }
    else {
        level=Level;
        if(previousLevels.size()==0 or level!=previousLevels.back())
            previousLevels.push_back(level);
    }
}

void PrintOutput::off(std::string Which){
    std::vector<std::string> allowed({"None","Both","Screen","File"});
    if(std::find(allowed.begin(),allowed.end(),Which)==allowed.end())
        DEVmessage("Illegal Which="+Which+", allowed:"+tools::str(allowed,0," "));

    if(Which=="File" or Which=="Both")out->remove("File");
    if(Which=="Screen" or Which=="Both")out->remove("Screen");
}

void PrintOutput::timerWrite(){
    PrintOutput::title("timer information");
    for(auto os: out->Streams())Timer::write(*os.second);
}

void PrintOutput::title(string Title){
    flush();
    *out<<'\n'<<'\n'<<" ";
    for(unsigned int k=0;k<3;k++)*out<<titleSeparator;
    if(Title.length()>0)*out<<"  "+Title+"  ";
    for(unsigned int k=Title.length()+8;k<lengthLine;k++)*out<<titleSeparator;
    *out<<'\n';
}
void PrintOutput::subTitle(string Sub){
    flush();
    *out<<" "<<Sub<<'\n';
}
void PrintOutput::verbatim(string Mess,ostream *Ostr){
    flush();
    if(level!="full")return;
    if(Ostr==0)*out<<Mess;
    else *Ostr<<Mess;
}

void PrintOutput::message(string Mess, ostream *Ostr, bool Final, int Limit){
    int cnt=messageCount[Mess];
    messageCount[Mess]++;
    if(Limit<0 and cnt>=-Limit)return;
    if(Limit>0 and cnt>= Limit)return;
    flush();
    if(level!="full")return;
    string line;
    if(Mess.length()<100)
        line=" *** "+Mess+" ***\n";
    else {
        line=Mess;
        size_t pos=0;
        while(string::npos!=(pos=line.find("\n",pos+10)))line.insert(pos+1,"   * ");
        line=" *** "+line+"\n";
    }

    // if message during progess monitor
    if(progressStatus!=""){
        line="\n"+line;
        progressStatus="";
    }

    if(Ostr==0)*out<<line;
    else *Ostr<<line;
    if(Final)finalMessage.push_back(line);
}
void PrintOutput::DEVmessage(std::string Mess, std::ostream *Ostr, bool Final){
#ifdef _DEVELOP_
    message("(Developer) "+Mess,Ostr,Final);
#endif
}
void PrintOutput::DEVwarning(std::string Mess, int Limit, std::ostream *Ostr){
#ifdef _DEVELOP_
    warning("(Developer) "+Mess,Limit,Ostr);
#endif

}

void PrintOutput::terminationMessages(ostream *Ostr){
    if(Ostr==0)for(int k=0;k<finalMessage.size();k++)*out<<"\n";
    else       for(int k=0;k<finalMessage.size();k++)*Ostr<<"\n";
    if(Ostr==0)for(int k=0;k<finalMessage.size();k++)*out<<finalMessage[k];
    else       for(int k=0;k<finalMessage.size();k++)*Ostr<<finalMessage[k];
}

void PrintOutput::matrix(const Eigen::MatrixXd & Mat,int Digits){

    PrintOutput::newRow();
    PrintOutput::rowItem("");
    double eps=std::pow(0.1,Digits);
    for(int l=0;l<Mat.cols();l++)PrintOutput::rowItem(int(l));
    for(int k=0;k<Mat.rows();k++){
        PrintOutput::newRow();
        PrintOutput::rowItem(tools::str(k)+"|");
        for(int l=0;l<Mat.cols();l++){
            PrintOutput::rowItem(std::abs(Mat(k,l))<eps?"  .  ":tools::str(Mat(k,l),Digits));
        }
    }
    PrintOutput::end();
}

void PrintOutput::progressStart(string Mess){
    progressHeader=Mess;
}
void PrintOutput::progressStop(string Mess, ostream *Ostr){
    if(progressStatus=="")return;
    if(Ostr==0)*out<<" "+Mess<<"\n";
    else      *Ostr<<" "+Mess<<"\n";
    progressStatus="";
}

void PrintOutput::progress(string Mess,ostream *Ostr){
    if(level=="off")return;
    string line=Mess;

    // first call after non-empty progressStart(...)
    if(progressHeader!="")line=progressHeader+": "+Mess;
    progressHeader="";

    progressStatus+=line;
    if(Ostr==0)*out<<line;
    else      *Ostr<<line;
    if(progressStatus.length()>80){
        if(Ostr==0)*out<<"\n";
        else      *Ostr<<"\n";
        progressStatus="";
    }
}

static std::map<std::string,int> _explanations;

void PrintOutput::warning(string Mess, int Limit, ostream *Ostr,std::string Explanations){
    int cnt=warningCount[Mess];
    warningCount[Mess]++;
    if(std::abs(cnt)==Limit){
        std::string line="                 !!! Limit="+tools::str(Limit)+" reached, no further warnings: "+Mess.substr(0,20)+"... !!!\n\n";
        if(Ostr==0)*out<<line;
        else *Ostr<<line;
    }

    if(Limit<0 and cnt>=-Limit)return;
    if(Limit>0 and cnt>= Limit)return;

    flush();
    // split into multi-line warning
    vector<string> mlin=tools::splitString(Mess,'\n');
    string line;
    for(size_t n=0;n<mlin.size();n++){
        line=mlin[n];
        if(n==0)line=" !!! WARNING --- "+line;
        else    line="             ... "+line;
        if(n==mlin.size()-1)line+=" --- WARNING !!!";
        line+="\n";
        if(Ostr==0)*out<<line;
        else *Ostr<<line;
    }
    if(Explanations!=""){
        _explanations[Mess.substr(0,40)]++;
        if(_explanations[Mess.substr(0,40)]==1){
            if(Ostr)*Ostr<<Explanations+"\n";
            *out<<Explanations+"\n\n";
        }
    }
    line="\n";
    if(Ostr==0)*out<<line;
    else *Ostr<<line;
}
void PrintOutput::warningList(){
    if(warningCount.size()==0)return;
    vector<string> list;
    list=tools::splitString(tools::listMapKeys(warningCount,"\f"),'\f');
    subTitle("\n !!! --- CURRENT WARNINGS --- !!!");
    for(unsigned int k=0;k<list.size();k++){
        std::string mess=list[k];
        std::replace(mess.begin(),mess.end(),'\n',' ');
        if(mess.length()>80)mess=mess.substr(0,80)+"...";
        lineItem(tools::str(warningCount[list[k]]),mess);
        newLine();
    }
    //    paragraph();
}

void PrintOutput::paragraph(){end();*out<<'\n';}
void PrintOutput::newLine(){end();*out<<'\n';}

void PrintOutput::lineItem(string Name, double Value, string Subtitle, unsigned int Digits){lineItem(Name,tools::str(Value,Digits,DBL_MAX/2),Subtitle);}
void PrintOutput::lineItem(string Name, int Value, string Subtitle)   {lineItem(Name,tools::str(Value),Subtitle);}
void PrintOutput::lineItem(string Name, string Value, string Subtitle){
    if(Subtitle.length()>0)subTitle(Subtitle);
    if(namesLine.length()!=0)
    {   if(Name.length()>0)namesLine +=nameSeparator;
        valuesLine+=itemSeparator;
    }
    namesLine +=Name;
    valuesLine+=Value;
}

//void PrintOutput::newRow(){row++;col=1;}
//void PrintOutput::rowItem(double Value,unsigned int Width){rowItem(tools::str(Value,Width,DBL_MAX/2));}
//void PrintOutput::rowItem(   int Value){rowItem(tools::str(Value));}
//void PrintOutput::rowItem(unsigned int Value){rowItem(tools::str(Value));}
//void PrintOutput::rowItem(string Value){
//    if(row==0 or col==0)ABORT("open newRow befor inserting rowItem");
//    for (size_t j=table.size();j<col;j++)table.push_back(vector<string>(0));
//    for (size_t i=table[col-1].size();i<row;i++)table[col-1].push_back("");
//    table[col-1][row-1]=Value;
//    col++;
//}

static std::shared_ptr<OutputTable> outputTable;
void PrintOutput::newRow(){
    if(not outputTable)outputTable.reset(new OutputTable());
    outputTable->newRow();
}
void PrintOutput::rowItem(double Value,unsigned int Width){outputTable->rowItem(Value,Width);}
void PrintOutput::rowItem(   int Value){outputTable->rowItem(Value);}
void PrintOutput::rowItem(unsigned int Value){outputTable->rowItem(size_t(Value));}
void PrintOutput::rowItem(string Value){outputTable->rowItem(Value);}

static int rowFlushed=-1;
void PrintOutput::end(){
    flush();
    outputTable.reset();
//    table.clear();
    rowFlushed=-1;
//    row=0;
//    col=0;
}

void PrintOutput::updateTable(){
    if(outputTable){
        auto lines=outputTable->format(rowFlushed+1);
        rowFlushed=outputTable->rows()-1;
        for(auto l: lines)*out<<l<<"\n";
    }
    else {
        rowFlushed=-1;
    }
}

//void PrintOutput::updateTable(){
//    // values in table
//    if(tableWidth.size()>0 or table.size()>0){
//        bool newTable=tableWidth.size()==0;
//        unsigned int length=0;
//        for (unsigned int n=0;n<table.size();n++){
//            if(newTable){
//                tableWidth.push_back(table[n][0].length());
//                for (unsigned int k=1;k<table[n].size();k++)
//                    tableWidth.back()=max(tableWidth.back(),(unsigned int)(table[n][k].length()));
//            }
//            length=max(length,(unsigned int)(table[n].size()));
//        }
//        for(unsigned int l=rowFlushed+1;l<length;l++){
//            for(unsigned int n=0;n<table.size();n++){
//                if(n==0)*out<<indent;
//                if(table[n].size()>l)*out<<setw(tableWidth[n])<<table[n][l];
//                else *out<<setw(tableWidth[n])<<" ";
//                if(n<table.size()-1)*out<<columnSeparator;
//            }
//            *out<<'\n';
//            if(newTable and l==0){
//                // first line in new table is header
//                *out<<indent;
//                for(unsigned int n=0;n<table.size();n++){
//                    if(n>0)*out<<" ";
//                    *out<<"  ";
//                    for(unsigned int k=1;k<tableWidth[n]-1;k++)*out<<"-";
//                }
//                *out<<'\n';
//            }
//        }
//        tableWidth.clear();
//        rowFlushed=length-1;
//        col=0;
//    }
//}

void PrintOutput::flush(){
    // values in a single line
    if(namesLine.length()>0){
        *out<<" "<<setw(max(namesTab,(unsigned int)(namesLine.length()+2)))<<namesLine+namesEnd+" ";
        *out<<valuesLine;
        namesLine="";
        valuesLine="";
    }

    // values in table
    updateTable();

}


