// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef STRINGAPPEND_H
#define STRINGAPPEND_H

#include <string>
#include <vector>
#include "stringTools.h"

#define SEP(separator) Str::sep(separator)
#define Sout Str(""," ")
#define Sstr Str(""," ")
#define Sendl Str::_Print(tools::str(__LINE__),std::string(__FILE__).substr(std::string(__FILE__).rfind("/")+1))

class Index;

/// \ingroup IO
/// @brief Convenience class for converting and appending to string
///
/// <br> Str's are concatenated by the + operator: Str(first,sep)+"other"+1 results in first+sep+"other"+sep+to_string(1)
/// <br> types are converted automatically upon append
/// <br> "+" is resolved from left to right, place brackets as needed for precedence:
/// e.g. Str::s+(7+8) gives "15", while Str::s+7+8 gives "7 8"
class Str : public std::string
{
    static std::string CurrentSeparator;

    int _wid;

public:
    ~Str(){}
    class _Print{
        std::string _pre;
        std::string _post;
    public:
        _Print():_pre(""),_post(""){}
        _Print(std::string Pre,std::string Post=""):_pre(Pre+": "),_post(Post){}
        std::string pre() const {return _pre;}
        std::string post() const {return _post;}
    };
    class File{
        friend class Str;
        std::string _name;
    public:
        File(std::string Name):_name(Name){}
    };
    class Separator{
        friend class Str;
        std::string _sep;
    public:
        Separator operator()(std::string Sep){_sep=Sep;return *this;}
    };

    static void Test();

    Str():_wid(0){}
    Str(std::string String,
        std::string Sep=CurrentSeparator/** how to separate upon concatenation */,
        int Width=0 /** width of enty, pad by blanks from the left, */);

    Str(char Char,
        std::string Sep=CurrentSeparator /** how to separate upon concatenation */,
        int Width=0 /** width of enty, pad by blanks from the left, */);

    /// append Arg as a string, use previous separator
    template<class T>  Str & operator+(T Arg){
        Str s(tools::str(Arg),CurrentSeparator,_wid);
        std::string ss=CurrentSeparator+s;
        *this=std::string::operator+=(ss);
        return *this;
    }
    Str & operator+(long double Arg){return operator+(double(Arg));}
    Str & operator+(std::string & Arg){return *this=std::string::operator+=(CurrentSeparator+Arg);}

    template<class T>
    Str & operator+=(T Arg){return operator+(Arg);}


    static Str emptyStr;  ///< (Macro STR) easy string for concatenating: Str::s+"something"+var+vec returns string with default separated items
    static _Print print; ///< (Macro PRINT) easy printing: Str("this string")+Str::print; prints to std::cout.
    static Separator sep; ///< (Macro SEP) change separator for further concatenations: Str("a","=")+x+Str::sep("...")+"y"+42 returns Str("a=x...y...42")
    void operator+(_Print Arg); ///< print to cout and clear, usage: Str::s,"my text",x,Str::print; or with macros: STR,"my text",x,PRINT;
    void operator+(File Arg);  ///< print to file and clear
    /// change separator, usage: Str("a","=")+x+Str::sep("...")+"y"+42 returns Str("a=x...y...42")
    Str & operator+(Separator Arg);

};
#endif // STRINGAPPEND_H
