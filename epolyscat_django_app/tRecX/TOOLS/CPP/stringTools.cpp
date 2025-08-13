// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "stringTools.h"
#include "tools.h"
#include "str.h"
#include "abort.h"
#include "algebra.h"
#include "Eigen/Dense"

namespace tools{

using namespace std;

size_t levenshtein_distance(const std::string S1, const std::string S2)
{
    size_t n=S1.size();
    size_t m=S2.size();
    ++n; ++m;
    size_t* d = new size_t[n * m];
    memset(d, 0, sizeof(size_t) * n * m);
    for (size_t i = 1, im = 0; i < m; ++i, ++im){
        for (size_t j = 1, jn = 0; j < n; ++j, ++jn){
            if (S1[jn] == S2[im]){
                d[(i * n) + j] = d[((i - 1) * n) + (j - 1)];
            }
            else{
                d[(i * n) + j] = min(d[(i - 1) * n + j] + 1, /* A deletion. */
                        min(d[i * n + (j - 1)] + 1, /* An insertion. */
                        d[(i - 1) * n + (j - 1)] + 1)); /* A substitution. */
            }
        }
    }
    size_t r = d[n * m - 1];
    delete [] d;
    return r;
}




/// base name of file
std::string fileBase(const std::string & File){
    size_t i0=File.rfind("/")+1;
    if(i0==std::string::npos)i0=0;
    return File.substr(i0,File.rfind(".")-i0);
}

/// remove leading blanks and returns from string
std::string lcropString(std::string s){while(s.length() and std::isspace(s.front()))s.erase(s.begin());return s;}
std::string rcropString(std::string s){while(s.length() and std::isspace(s.back() ))s.pop_back();return s;}
std::string cropString(std::string s){return lcropString(rcropString(s));}

string toLower(const string Inp) {
    string out(Inp);
    for (int b=0; b < Inp.length(); b++)out[b]=tolower(Inp[b]);
    return out;
}

std::string unquote(string s){
    // single quotes
    size_t i0=s.find("'"),i1=s.rfind("'");
    if(i0==0 and i1==s.length()-1)return s.substr(1,s.length()-2);
    i0=s.find('"');
    i1=s.rfind('"');
    if(i0==0 and i1==s.length()-1)return s.substr(1,s.length()-2);
    return s;
}

std::string str(std::string S){return S;}

std::string str(bool Bool){
    if(Bool)return "true";
    return "false";
}

/// convert Number to std::string (overloaded)
string str(char Arg, int Len){
    ostringstream oss;
    if(Len>0){oss<<std::setw(Len)<< Arg;return oss.str();}
    else     {oss<<                 Arg;return oss.str();}
}
string str(string Arg, int Len){
    ostringstream oss;
    if(Len>0){oss<<std::setw(Len)<< Arg;return oss.str();}
    else     {oss<<                 Arg;return oss.str();}
}
string str(const int          Number, int Len, char Fill){
    ostringstream oss;
    if(Len>0){oss<<std::setw(Len)<<setfill(Fill)<< Number;return oss.str();}
    else     {oss<<                                Number;return oss.str();}
}
string str(const unsigned int Number, int Len, char Fill){
    ostringstream oss;
    if(Len>0){oss<<std::setw(Len)<<setfill(Fill)<< Number;return oss.str();}
    else     {oss<<                                Number;return oss.str();}
}
string str(const size_t Number, int Len, char Fill){
    ostringstream oss;
    if(Len>0){oss<<std::setw(Len)<<setfill(Fill)<< Number;return oss.str();}
    else     {oss<<                                Number;return oss.str();}
}
string str(const long int Number, int Len, char Fill){
    ostringstream oss;
    if(Len>0){oss<<std::setw(Len)<<setfill(Fill)<< Number;return oss.str();}
    else     {oss<<                                Number;return oss.str();}
}
std::string str(std::complex<double>* Pointer, int Dum){
    ostringstream oss;oss<<Pointer;return oss.str();
}

string str(const size_t Size){ 
    ostringstream ss;
    ss<<Size;
    return ss.str();
}
string str(const long int Size){
    ostringstream ss;
    ss<<Size;
    return ss.str();
}
string str(double Number, int Prec, double Inf){
    if(Prec==0)return string(1,tools::zero(complex<double>(Number)));
    if(Inf>0.){
        string inf=Algebra::Infty;
        inf=inf.substr(0,Prec);
        if(Number>= Inf)return     inf;
        if(Number<=-Inf)return "-"+inf;
    }
    ostringstream oss; oss<<std::setprecision(Prec)<< Number;return oss.str();}
string str(complex<double> Number, int Prec){if(Prec==0)return string(1,tools::zero(Number));return "("+str(real(Number),Prec)+","+str(imag(Number),Prec)+")";}

std::string str(const Eigen::MatrixXcd & M, int Digits){
    ostringstream oss;
    // determine maximal width
    int maxW=1;
    if(Digits>0)
        for(int j=0;j<M.cols();j++)
            for(int i=0;i<M.rows();i++)
                maxW=std::max(maxW,(int)tools::str(M(i,j),Digits).length()+1);

    for(int i=0;i<M.rows();i++){
        for(int j=0;j<M.cols();j++)
            oss<<std::setw(maxW)<<tools::str(M(i,j),Digits);
        oss<<std::endl;
    }

    return oss.str();
}

/// number of (non-overlapping) occurrences of substring in String
unsigned int subStringCount(const std::string String,const std::string Sub){
    if(Sub.length()==0)return 0;
    unsigned int count=0;
    size_t pos=String.find(Sub,0);
    while(pos!=std::string::npos){
        count++;
        pos=String.find(Sub,pos+Sub.length());
    }
    return count;
}
/// number of (non-overlapping) occurrences of substring in String
unsigned int subStringCount(const std::string String,const std::string Sub, string Left, string Right){
    if(Sub.length()==0)return 0;
    unsigned int count=0;
    size_t pos=findFirstOutsideBrackets(String,Sub,Left,Right,0);
    while(pos!=std::string::npos){
        count++;
        pos=findFirstOutsideBrackets(String,Sub,Left,Right,pos+Sub.length());
    }
    return count;
}
/// replace in String first appearance of OldSub with NewSub
std::string substringReplace(const std::string & S, const std::string & OldSub, const std::string & NewSub){
    size_t pos=S.find(OldSub);
    std::string s(S);
    if(pos!=std::string::npos){
        s=S.substr(0,pos)+NewSub+S.substr(pos+OldSub.length());
    }
    return s;
}
/// replace in String all appearances of OldSub with NewSub
std::string substringReplaceAll(const std::string & S, const std::string & OldSub, const std::string & NewSub){
    std::string s(S);
    while(s.find(OldSub)!=std::string::npos)
        s=substringReplace(s,OldSub,NewSub);
    return s;
}

/// convert strings to numbers
int    string_to_int(const string &Text){ return atoi(Text.c_str()); }

double string_to_double(const string &Text,const double Inf){
    if(tools::cropString(Text)==Algebra::Infty)return  Inf;
    if(tools::cropString(Text)=="-"+Algebra::Infty)return -Inf;
    char* pEnd;
    return strtod(Text.c_str(),&pEnd); }

complex<double> string_to_complex(const string &Text){
    char* pEnd;
    size_t comma=Text.find(",");
    if(comma==string::npos)return complex<double>(strtod(Text.c_str(),&pEnd),0);
    else return complex<double>(strtod(Text.substr(0,comma).c_str(),&pEnd),strtod(Text.substr(comma+1).c_str(),&pEnd));
}

//// conversion to bool: must be either true or false
bool string_to_bool(const string &Text){
    if (Text != "true"){
        if (Text != "false")ABORT("bool string must be true or false, is: " + Text);
        return false;
    }
    return true;
}


/// [beg,end,pts] ... expands to { beg, beg+(end-beg)/(pts-1),..., end} <br>
/// [beg] or beg or beg==end ... single point
//std::vector<double> rangeToGrid(string Range, int Points){
//    vector<string> parts=tools::splitString(tools::stringInBetween(Range,"[","]"),',');
//    vector<double> def;
//    for(string p: parts)def.push_back(Algebra::realConstant(p));

//    if(def.size()==0)return {0.,0.,0};
//    if(def.size()==1 or def[0]==def[1])return {def[0]};

//    int points=Points;
//    if(def.size()==3)points=max(2,int(def[2]));

//    vector<double>grid;
//    for(int k=0;k<points;k++)grid.push_back(def[0]+k*(def[1]-def[0])/double(points-1));
//    return grid;
//}
/// get equidistant grid from Range string (should replace the stringTools version), attach units like' eV' or '~eV'
std::vector<double> rangeToGrid(std::string Range /** [beg,end,pts] */,
                                int Points        /** default for pts */,
                                std::vector<std::string> Delim /** alternate deliminators for Range */)
{
    if(Delim.size()!=3)DEVABORT("need format Delim={begString,sepChar,endString}");
    if(Delim[0]!="" and Delim[1]!="")Range=tools::stringInBetween(Range,"[","]");
    if(Delim[1].length()!=1)DEVABORT(Sstr+"can only use single-character separator, found: Delim="+Delim);
    vector<string> parts=tools::splitString(Range,Delim[1][0]);
    vector<double> def;
    for(string p: parts)def.push_back(Algebra::realConstant(p));

    if(def.size()==0)return {0.,0.,0};
    if(def.size()==1 and def[0]==def[1])return {def[0]};

    int points=Points;
    if(def.size()==3)points=max(2,int(def[2]));

    vector<double>grid;
    for(int k=0;k<points;k++)grid.push_back(def[0]+k*(def[1]-def[0])/double(points-1));
    return grid;
}


std::vector<int> string_to_intVec(std::string S){
    std::vector<int> ir;
    std::string s(cropString(S));
    if(s.find(':')==std::string::npos){
        if(s[0]=='{')s=s.substr(1);
        if(s.back()=='}')s.pop_back();
        if(s.find('{')!=std::string::npos)goto Terminate;
        if(s.find('}')!=std::string::npos)goto Terminate;
        std::vector<std::string> sv(tools::splitString(s,','));
        for (auto v: sv)ir.push_back(tools::string_to_int(v));
    }
    else {
        std::vector<std::string> sv(tools::splitString(s,':'));
        int  ibeg=tools::string_to_int(sv[0]);
        int  iend=tools::string_to_int(sv[1]);
        int  inc=sv.size()<3?1:tools::string_to_int(sv[2]);
        if(inc==0)ABORT("cannot have 0 increment for integer range, got: "+S);
        for(int k=ibeg;(inc>0 and k<iend) or (inc<0 and k>iend);k=k+inc)ir.push_back(k);
    }
    if(ir.size()==0)goto Terminate;
    return ir;
Terminate:
    ABORT("Cannot inteprete integer range "+S+"\n admissable format(s): 2,1,6,7,.. or {2,1,6,7,..}  or i0:i1 or i0:i1:inc");
    return {};
}

std::vector<double> string_to_doubleVec(std::string S){
    std::vector<double> ir;
    std::string s(cropString(S));
    if(S.front()=='{' and S.back()=='}'){
        s.pop_back();
        s=s.substr(1);
        if(s.find(':')!=std::string::npos)goto Terminate;
        if(s.find('{')!=std::string::npos)goto Terminate;
        if(s.find('}')!=std::string::npos)goto Terminate;
        std::vector<std::string> sv(tools::splitString(s,','));
        for (auto v: sv)ir.push_back(tools::string_to_double(v));
    }
    else if(s.find(':')!=std::string::npos){
        std::vector<std::string> sv(tools::splitString(s,':'));
        if(sv.size()<2)
            ir.push_back(tools::string_to_int(sv[0]));
        else {
            int  ibeg=tools::string_to_int(sv[0]);
            int  iend=tools::string_to_int(sv[1]);
            int  inc=sv.size()<3?1:tools::string_to_int(sv[2]);
            if(inc==0)ABORT("cannot have 0 increment for integer range, got: "+S);

            for(int k=ibeg;(inc>0 and k<iend) or (inc<0 and k>iend);k=k+inc)ir.push_back(double(k));
        }
    }
    return ir;
Terminate:
    ABORT("Cannot inteprete integer range "+S+"\n admissable format(s): {2,1,6,7,..}, i0:i1, i0:i1:inc");
    return {};
}

std::vector<std::string> splitString(const std::string s, char delimiter, const string LeftBracket, const string RightBracket) {
    string s0(s);
    // remove white-space, if blank delimiter
    if(delimiter==' ')s0=tools::cropString(s);
    if(s0=="")return {};

    std::vector<std::string> elements,sep;
    splitString(s0,string(1,delimiter),elements,sep,LeftBracket,RightBracket);

    // remove white-space, if blank delimiter
    if(delimiter==' '){
        for(string &e: elements)e=tools::cropString(e);
        if(elements.size() and elements.back()=="")elements.pop_back();
    }

    if(elements.size()==0 and s!="")elements.push_back(s);
    return elements;
}

std::string abbreviate(const std::string &LongString, size_t Length){
    if(LongString.length()<=Length)return LongString;
    if(LongString.length()<=Length+5)return LongString.substr(0,Length-3)+"...";
    return LongString.substr(0,Length/2-2)
            +" ... "
            +LongString.substr(LongString.length()-Length/2+3);
}

std::vector<std::string> splitString(const std::string &String, const std::string &Deliminator){
    std::vector<std::string> res;
    size_t beg=0;
    for(size_t end=String.find(Deliminator,beg);end!=std::string::npos;end=String.find(Deliminator,end+1)){
        res.push_back(String.substr(beg,end-beg));
        beg=end+Deliminator.length();
    }
    res.push_back(String.substr(beg));
    return res;
}

std::string joinString(const std::vector<std::string> & VecStrings, const string &Deliminator){
    std::string res;
    for(auto s: VecStrings)res+=s+Deliminator;
    return res.substr(0,res.length()-Deliminator.length());
}


/// next occurrence of any character of Sub that is not between any pair of the brackets
size_t OLDfindOutsideBrackets(bool First, const string S, const string &Sub, const string &Left, const string &Right, size_t Pos){
    if(Left.length()!=Right.length()){cout<<"ERROR: number of left brackets does not match right\n";exit(1);}
    size_t i0;
    if(First)i0=S.find_first_of(Sub,Pos);
    else     i0=S.find_last_of(Sub,Pos);
    if(i0==string::npos)return i0;
    int openBracket=0;
    for (unsigned int k=0;k<Left.length();k++){
        string L=Left.substr(k,1),R=Right.substr(k,1);
        // it may happen that delimiters at one side occur multiple times, count that
        unsigned int Lmult=subStringCount(Left,L),Rmult=subStringCount(Right,R);
        if(L!=R)openBracket+=subStringCount(S.substr(0,i0),L)*Rmult-subStringCount(S.substr(0,i0),R)*Lmult;
        else openBracket+=subStringCount(S.substr(0,i0),L)%2;
    }
    if(openBracket!=0){
        if(i0+1==Pos)return string::npos; // no position outside brackets
        return findOutsideBrackets(First,S,Sub,Left,Right,i0+1);
    }
    return i0;
}
size_t findOutsideBrackets(bool First, const string S, const string &Sub, const string &Left, const string &Right, size_t Pos){
    if(Left.length()!=Right.length()){cout<<"ERROR: number of left brackets does not match right\n";exit(1);}
    size_t i0;
    if(First)i0=S.find_first_of(Sub,Pos);
    else     i0=S.find_last_of(Sub,Pos);
    if(i0==string::npos)return i0;

    std::string s(S);
    int openBracket=0;
    int closeBracket=0;
    for (unsigned int k=0;k<Left.length();k++){
        string L=Left.substr(k,1),R=Right.substr(k,1);

        // it may happen that delimiters at one side occur multiple times, count that
        unsigned int Lmult=subStringCount(Left,L),Rmult=subStringCount(Right,R);
        if(L!=R)openBracket+=subStringCount(S.substr(0,i0),L)*Rmult-subStringCount(S.substr(0,i0),R)*Lmult;
        else    openBracket+=subStringCount(S.substr(0,i0),L)%2;
        if(L!=R)closeBracket+=subStringCount(S.substr(i0+1),R)*Rmult-subStringCount(S.substr(i0+1),L)*Lmult;
        else    closeBracket+=subStringCount(S.substr(i0+1),R)%2;
        closeBracket=max(0,closeBracket);
    }
    if(openBracket!=0 and closeBracket!=0){
        if(i0+1==Pos and openBracket==closeBracket)return string::npos; // not position outside brackets
        return findOutsideBrackets(First,S,Sub,Left,Right,i0+1);
    }
    return i0;
}


size_t findFirstOutsideBrackets(const string &S,const string Sub,const string &Left, const string &Right, size_t Pos){return findOutsideBrackets(true,S,Sub,Left,Right,Pos);}
size_t findLastOutsideBrackets (const string &S,const string Sub,const string &Left, const string &Right, size_t Pos){return findOutsideBrackets(false,S,Sub,Left,Right,Pos);}

/// split a string into substrings at a set of separators, return actual separators in vector
/// first item can, but does not need to be preceded by a separator
/// e.g. split  s="a + b - ccc " with Separators="+-" into Elem=("a "," b "," ccc "), Sep=(" ","+","-")
/// any separators inside a LeftBracket/RightBracket pair will be ignored
/// this allows splitting expressions <dJd+qJq>+<J> into <dJd+qJq> and <J>
void splitString(const string S, string Separators, vector<string> &Elem, vector<string> & Sep,const string LeftBracket,const string RightBracket) {
    Elem.clear();
    Sep.clear();
    string s;

    // if blank separation is an option, reduce to single blank per interval
    if(Separators.find(" ")!=string::npos){
        unsigned int k0=0,k1=S.find(" ");
        while(k0<S.length()){
            s+=tools::cropString(S.substr(k0,k1-k0))+" ";
            k0=S.find_first_not_of(" ",k1);
            k1=S.find(" ",k0);
            if(k1>S.length())k1=S.length();
        }
        tools::cropString(s);
    } else {
        s=tools::cropString(S);
    }
    size_t i0=0,i1=findFirstOutsideBrackets(s,Separators,LeftBracket,RightBracket);
    Elem.push_back(cropString(s.substr(i0,i1)));

    if(Elem.back()=="")Elem.pop_back(); // nothing before first separator
    else Sep.push_back(" ");
    while(i1!=std::string::npos){
        i0=i1;
        i1=findFirstOutsideBrackets(s,Separators,LeftBracket,RightBracket,i1+1);
        Elem.push_back(s.substr(i0+1,i1-i0-1));
        Sep.push_back(s.substr(i0,1));
    }
    if(Elem.size()==0 and S!="")Elem.push_back(S);
}

/// Left and Right constain subsequent pairs of brackets, e.g. ('| and )'> for (...), '...', |...>
/// if no matching pair is found, the original string is returned
string stringInBetween(string S, string Left, string Right,bool Full){
    if(Full){
        size_t il=S.find(Left),ir=S.rfind(Right);
        if(il<ir and ir!=string::npos)return S.substr(il+Left.length(),ir-il-Left.length());
    } else {
        for(unsigned int k=0;k<min(Left.size(),Right.size());k++){
            size_t il=S.find(Left.substr(k,1)),ir=S.rfind(Right.substr(k,1));
            if(il<ir and ir!=string::npos)return S.substr(il+1,ir-il-1);
        }
    }
    return S;
}

std::string stringStrip(std::string Str, std::string End){
    return Str.find(End)+End.length()==Str.length()?Str.substr(0,Str.find(End)):Str;
}

}




