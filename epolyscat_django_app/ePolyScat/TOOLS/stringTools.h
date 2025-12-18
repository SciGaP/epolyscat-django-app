// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef STRINGTOOLS_H
#define STRINGTOOLS_H

#include "toolsHeader.h"
#include "qtEigenDense.h"
#include <set>

namespace tools{

/** @defgroup  StringTools
 *  @ingroup Tools
 *  \brief String-specific tools
 * @{
 */

using namespace std;

/// base name of file
std::string fileBase(const std::string & File);


std::string lcropString(std::string s);///< remove leading blanks from string
std::string rcropString(std::string s);///< remove trailing blanks from string
std::string cropString(std::string s);///< return without white outer white space
std::string unquote(std::string s); ///< return without pair of quotes  at first and last position in string

//void splitString(const std::string s,char delimiter, std::vector<std::string>& elements);
//std::vector<std::string> splitString(const std::string s, char delimiter);

/// convert Number to std::string (overloaded)
std::string str(char Arg, int Len=0);
std::string str(std::string Arg, int Len=0);
std::string str(double Number, int Prec=14, double Inf=0.);
std::string str(std::complex<double> Number,int Prec=3);
std::string str(const int          Number, int Len=0, char Fill=' ');
std::string str(const unsigned int Number, int Len=0, char Fill=' ');
std::string str(const size_t Size, int Len=0, char Fill=' ');
std::string str(const long int Size, int Len=0, char Fill=' ');
std::string str(std::complex<double>* Pointer, int dum=0);
std::string str(bool Bool);
std::string str(const Eigen::MatrixXcd & M, int Digits=1);//{ostringstream oss; oss<< M;return oss.str();}

std::string toLower(const std::string Inp);

template<typename T>std::string str(T*P, int dum=0){ostringstream oss; oss<< P;return oss.str();}
//template<typename T>std::string str(const T&P, int dum=0){ostringstream oss; oss<< P;return oss.str();}

template<typename T> std::string str(const std::vector<T> & vec, int Prec=3, string Sep=" "){
    string s;
    for(unsigned int n=0;n<vec.size();n++){if(n==0)s=str(vec[n],Prec); else s+=Sep+str(vec[n],Prec);};
    return s;
}

template<typename T> std::string str(const std::set<T> & vec, int Prec=3, string Sep=" "){
    string s;
    for(auto p=vec.begin();p!=vec.end();p++){if(p==vec.begin())s=str(*p,Prec); else s+=Sep+str(*p,Prec);};
    return s;
}

/// get equidistant grid from Range string
std::vector<double> rangeToGrid(std::string Range /** [beg,end,pts] */,
                                int Points=11 /** default for pts */,
                                std::vector<std::string> Delim={"[",",","]"} /** alternate deliminators for Range */ );

/// @brief vector of all keys in a map [works for keys, where str(class K) is defined]
template <class K, class V>
std::vector<std::string> vectorMapKeys(std::map<K,V> Map){
    std::vector<std::string> vKey;
    for(typename std::map<K,V>::iterator it=Map.begin();it!=Map.end();it++)vKey.push_back(str(it->first));
    return vKey;
}

template <class K, class V>
std::string listMapKeys(std::map<K,V> Map,
                        const std::string Sep=" ",/**< separator for the key strings */
                        const K* Except=0 /**< point to Key that should not be listed */
        ){
    std::string s;
    for(typename std::map<K,V>::iterator it=Map.begin();it!=Map.end();it++)
        if(Except==0 or *Except!=it->first)s+=str(it->first)+Sep;
    return s.substr(0,s.rfind(Sep));
}
/// @brief list of all keys in a map [works for keys, where str(class K) is defined]
template <class K, class V>
std::string listMap(std::map<K,V> Map,const std::string Sep="\n"/**< separator for strings and values */)
{
    std::string s;
    for(typename std::map<K,V>::iterator it=Map.begin();it!=Map.end();it++)
        s+=str(it->first)+"="+str(it->second)+Sep;
    return s.substr(0,s.rfind(Sep));
}

/// @brief Levenshtein distance - measure for similarity of strings
size_t levenshtein_distance(const string S1, const string S2);


/// @brief check for existence of Key in Map
template <class K, class V>
bool hasKey(std::map<K,V> Map,const K Key){
    for(typename std::map<K,V>::iterator it=Map.begin();it!=Map.end();it++)
        if(it->first==Key)return true;
    return false;
}

/// number of (non-overlapping) occurances of substring in String
unsigned int subStringCount(const std::string String,const std::string Sub);
unsigned int subStringCount(const std::string String,const std::string Sub,std::string Left,std::string Right);

/// replace in String first appearance of OldSub with NewSub (std::regex_replace did not function on Intel)
std::string substringReplace(const std::string & S, const std::string & OldSub, const std::string & NewSub);
/// replace in String all appearances of OldSub with NewSub (std::regex_replace did not function on Intel)
std::string substringReplaceAll(const std::string & S, const std::string & OldSub, const std::string & NewSub);


/// convert strings to numbers
int    string_to_int(const string &Text);
double string_to_double(const string &Text,const double Inf=DBL_MAX);
complex<double> string_to_complex(const string &Text);
/// conversion to bool: must be either true or false
bool string_to_bool(const string &Text);

/// strings to vectors as "{1,2,4,5}"->{1,2,4,5}, "3:5"->{3,4}, "9:3:-2"->{9,7,5}
std::vector<int> string_to_intVec(std::string S);

/// strings to vector (see string_to_intVec)
std::vector<double> string_to_doubleVec(std::string S);

/// next occurrence of any character of Sub that is NOT between any pair of the brackets
size_t findOutsideBrackets(bool First, const string S, const string &Sub, const string &Left, const string &Right, size_t Pos=0);
size_t findFirstOutsideBrackets(const string &S,const string Sub,const string &Left, const string &Right, size_t Pos=0);
size_t findLastOutsideBrackets (const string &S,const string Sub,const string &Left, const string &Right, size_t Pos=std::string::npos);

/// return string from after first occurrance of left until last occurance of right
/// Full true: match full strings Left/Right
///      false: match each pair of characters in Left/Right
std::string stringInBetween(std::string S, std::string Left, std::string Right, bool Full=false);

/// remove End from Str, return Str if it does not end in End
std::string stringStrip(std::string Str, std::string End);

/// split a string into substrings at a set of separators, return actual separators in vector
/// first item can, but does not need to be preceded by a separator
/// e.g. split  s="a + b - ccc " with Separators="+-" into Elem=("a "," b "," ccc "), Sep=(" ","+","-")
/// any separators inside a LeftBracket/RightBracket pair will be ignored
/// this allows splitting expressions <dJd+qJq>+<J> into <dJd+qJq> and <J>
void splitString(const string s, string Separators, vector<string> &Elem, vector<string> & Sep,
                 const string LeftBracket="",const string RightBracket="");

/// split by character, except inside bracket-pairs; delimiter=blank (' ') is special: remove any extra white-space
std::vector<std::string> splitString(const std::string s, char delimiter,
                                     const string LeftBracket="",const string RightBracket="");

/// if LongString exceeds Length, return beginning and end with ... in between
std::string abbreviate(const std::string &LongString, size_t Length=80);

/// split at Deliminator string, omit Delimitnator section
std::vector<std::string> splitString(const std::string &String, const std::string &Delimiter);
/// join into string using Deliminator between parts
std::string joinString(const std::vector<std::string> & VecStrings, const std::string &Delimiter);
/** @} */

}
#endif // STRINGTOOLS_H
