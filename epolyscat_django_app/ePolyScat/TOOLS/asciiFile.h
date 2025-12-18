// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ASCIIFILE_H
#define ASCIIFILE_H


#include <string>
#include <vector>
#include "limits.h"

#include <limits>       // std::numeric_limits
#include <stdio.h>      /* fopen, fputs, fclose, stderr */
#include <iostream>
#include <fstream>
#include <sstream>
#include "abort.h"
#include "folder.h"
#include "stringTools.h"
#include "qtEigenDense.h"

/// \ingroup IO
/// \brief open, write columns, rows, copy etc.
class AsciiFile
{
    static std::string comment;
    std::string _name;
    unsigned int _blank;
    std::string _sep;
public:
    AsciiFile(){}
    AsciiFile(std::string Name,
              unsigned int BlankInsert=INT_MAX /** insert blank if value at column number changes (gnuplot style) */,
              std::string Separator=", ",
              bool FixedLength=true);
    std::string name() const {return _name;}
    bool empty() const; ///< true if nothing in file

    void blankIfChange(int Col){_blank=Col;} ///< insert blank line when this Col changes value

    /// add a columns to existing file, requirse match of indicated two columns
    bool addCols(const std::vector<std::vector<double> > &Cols, std::vector<size_t> Match={0,0});

    /// read columns
    void readCols(std::vector<std::vector<double> > &Cols,
                  std::vector<unsigned int> NCols/** list of columns to read (empty=all)*/,
                  const std::string Sep /** column separator */,
                  std::vector<int> & DataSize /** size of each data set as found on file */,
                  const std::vector<int> DataSet=std::vector<int>()/** list of blank-line separated data sets to read (empty == all) */);

    /// read comma-separated or whitespace separated columns of numbers into Cols, NCols allows specifing column numbers, default is comma-separated
    void readCols(std::vector<std::vector<double> > &Cols,
                  std::vector<unsigned int> NCols=std::vector<unsigned int>()/** list of columns to read (default=all)*/,
                  const std::string Sep="," /** column separator */,
                  int DataSet=0 /** read DataSet'th blank separated set, ==-1 read all */);

    /// read multiple data sets (blank-row separated) into matrices
    void readCols(std::vector<Eigen::MatrixXd> &Cols,
                  std::vector<unsigned int> NCols=std::vector<unsigned int>()/** list of columns to read (default=all)*/,
                  const std::string Sep="," /** column separator */,
                  bool SetsRowWise=true /** each set will be one row of matrix */);

    /// read at single column
    void readCol(std::vector<double> &Col, unsigned int Pos=0 /** read Pos'th column */);

    void readComments(std::vector<std::string> &Comment);
    void writeComments(const std::vector<std::string> &Comment);
    /// insert comments into existing file
    void insertComments(const std::vector<std::string> &Comment,int Line=-2 /** Line<0: at line len(comments)+Line*/);

    void writeBlankRow(int Count=1) const;

    void writeCols(const std::vector<std::vector<std::complex<double>>> &Cols, unsigned int Digits=7) const;
    void writeCols(const std::vector<std::vector<double> > &Cols, unsigned int Digits=7) const;
    void writeRow(const std::vector<double> &Row, unsigned int Digits=7) const;
    void copy(std::string NewFile, bool Overwrite=false);
    void copy(std::ostream & Stream);

    /// acscii writes, mostly for trivially portable filing
    template<class T> void writeTagged(std::string Tag,const T & Data); ///< add line as Tag<<": "<<Data
    template<class T> void writeTagged(std::string Tag,const std::vector<T> & Data); ///< add line as Tag<<": "<<Data
    template<class T> void write( T * Data,size_t Size);

    /// the matching reads
    template<class T> void readTagged(std::string Tag, T & Data);
    template<class T> void readTagged(std::string Tag, std::vector<T> & Data);
    template<class T> void read( T * Data,size_t Size);
};

template<class T> void AsciiFile::writeTagged(std::string Tag,const T & Data){
    std::ofstream stream;
    stream.open(_name.c_str(),std::ios_base::out|std::ios_base::app);
    if(not stream.is_open())ABORT("could not open output file '"+_name+"'");

    stream << std::setprecision( std::numeric_limits<T>::digits10+2);
    stream<<Tag<<": "<<Data<<std::endl;
    stream.close();
}
template<class T> void AsciiFile::readTagged(std::string Tag,T & Data){
    if(not folder::exists(_name))ABORT("input file '"+_name+"' does not exist");
    std::ifstream stream;
    stream.open(_name.c_str(),std::ios_base::in);
    if(not stream.is_open())ABORT("could not open input file '"+_name+"'");

    // not very elegant...
    std::string tag;
    std::streampos pos;
    while(std::getline(stream,tag) and tag.find(Tag)!=0)pos=stream.tellg();
    stream.seekg(pos);

    stream>>tag>>Data;
    if(stream.fail())ABORT("error reading tag \""+Tag+"\" from '"+_name+"'");
    stream.close();
}

template<class T> void AsciiFile::writeTagged(std::string Tag,const std::vector<T> & Data){
    std::ofstream stream;
    stream.open(_name.c_str(),std::ios_base::out|std::ios_base::app);
    if(not stream.is_open())ABORT("could not open output file '"+_name+"'");

    stream << std::setprecision( std::numeric_limits<T>::digits10+2);
    stream<<Tag<<"["<<Data.size()<<"]:";
    for(int k=0;k<Data.size();k++)stream<<" "<<Data[k];
    stream<<std::endl;
    stream.close();
}

template<class T> void AsciiFile::readTagged(std::string Tag,std::vector<T> & Data){
    if(not folder::exists(_name))ABORT("input file '"+_name+"' does not exist");
    std::ifstream stream;
    stream.open(_name.c_str(),std::ios_base::in);
    if(not stream.is_open())ABORT("could not open input file '"+_name+"'");

    // not very elegant...
    std::string tag;
    std::streampos pos;
    while(std::getline(stream,tag) and tag.find(Tag)!=0)pos=stream.tellg();
    stream.seekg(pos);

    stream>>tag;
    Data.resize(tools::string_to_int(tools::stringInBetween(tag,"[","]")));
    for(int k=0;k<Data.size();k++)stream>>Data[k];
    if(stream.fail())ABORT("error reading tag \""+Tag+"\" from '"+_name+"'");
    stream.close();
}

template<class T> void AsciiFile::write( T * Data,size_t Size){
    std::ofstream stream(_name.c_str(),std::ios::out);
    if(not stream.is_open())ABORT("could not open output file '"+_name+"'");

    stream<<std::setprecision(17);
    for(size_t k=0;k<Size;k++)stream<<Data[k]<<" ";
    stream<<std::endl<<std::flush;
    stream.close();
}

template<class T> void AsciiFile::read( T * Data,size_t Size){
    if(not folder::exists(_name))ABORT("input file '"+_name+"' does not exist");
    std::ifstream stream(_name.c_str(),std::ios_base::in);
    if(not stream.is_open())ABORT("could not open input file '"+_name+"'");

    for(size_t k=0;k<Size;k++)stream>>Data[k];
    if(stream.fail())ABORT("read error on '"+_name+"'");
    stream.close();
}
#endif // ASCIIFILE_H
