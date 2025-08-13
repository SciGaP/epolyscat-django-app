// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "str.h"
#include "asciiFile.h"
#include "useMatrix.h"
#include "mpiWrapper.h"
#include "debugInfo.h"

using namespace std;
string Str::CurrentSeparator=" ";
Str::_Print Str::print;
Str::Separator Str::sep;
Str Str::emptyStr(""," ");

Str::Str(std::string String, std::string Sep, int Width):_wid(Width){
    CurrentSeparator=Sep;
    if(_wid!=0 and _wid>String.length())String=string(' ',_wid-String.length())+String;
    assign(String);
}

Str::Str(char Char, std::string Sep, int Width):_wid(Width){
    CurrentSeparator=Sep;
    string s(1,Char);
    if(_wid!=0)s=string(' ',_wid-1)+s;
    assign(s);
}

void Str::operator+(_Print Arg) {
    if(debug_tools::verboseLevel()<2)return;
    std::string mess=Arg.pre()+*this+" (Str::print) "+Arg.post();
    if(MPIwrapper::Size(MPIwrapper::worldCommunicator())==1)
        std::cout<<mess<<std::endl<<std::flush;
    else
        std::cout<<"<"<<MPIwrapper::Rank(MPIwrapper::worldCommunicator())<<"> "<<mess<<std::endl<<std::flush;

    clear();
    CurrentSeparator=" "; // restore default
}

Str & Str::operator+(Separator Arg) {
    CurrentSeparator=Arg._sep;
    return *this;
}

void Str::operator+(File F){
    ofstream stream;
    stream.open(F._name.c_str(),std::ios_base::out|std::ios_base::app);
    if(not stream.is_open())ABORT("could not open input file '"+F._name+"'");
    stream<<*this<<endl;
    clear();
    CurrentSeparator=" "; // restore default
}

void Str::Test(){
    complex<double>x(1.,1.);
    vector<int> i(3,6);
    Str("a","=")+x+Str::sep("...")+42+7+i+Str::print;
    Sstr+"try"+SEP("----")+"this"+Sendl;
    Sstr+"show a pointer"+SEP(" &")+&i+Sendl;
}
