// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "asciiFile.h"
#include "abort.h"
#include "folder.h"

#include "tools.h"
#include "str.h"

using namespace std;

string AsciiFile::comment="/#";

AsciiFile::AsciiFile(std::string Name, unsigned int BlankInsert, string Separator, bool FixedLength)
    :_name(Name),_blank(BlankInsert),_sep(Separator){if(not FixedLength)_sep=", ";}

void AsciiFile::readComments(std::vector<string> &Comment){
    ifstream stream(_name.c_str(),std::ios::in);
    string line;
    Comment.clear();
    while(getline(stream,line))
        if(line.find_first_of(comment)==0)
            Comment.push_back(line.substr(1));
}
void AsciiFile::writeComments(const std::vector<string> &Comment){
    ofstream stream(_name.c_str(),std::ios::out);
    for(unsigned int k=0;k<Comment.size();k++){
        if(Comment[k].find("#")!=0)stream<<"# ";
        stream<<Comment[k]<<endl;
    }
    stream<<flush;
}

void AsciiFile::insertComments(const std::vector<string> &Comment, int Line){
    std::vector<std::string> com;
    std::vector<std::vector<double>>cols;
    readComments(com);
    readCols(cols);
    if(com.size()==0)
        com=Comment;
    else{
        if(com.size()<size_t(abs(Line)))Line=com.size()-1;
        else if (Line<0)Line=com.size()+Line;
        com.insert(com.begin()+Line,Comment.begin(),Comment.end());
    }
    writeComments(com);
    writeCols(cols);
}

void AsciiFile::writeCols(const std::vector<std::vector<std::complex<double>>>& Cols, unsigned int Digits) const {
    std::vector<std::vector<double>>cDouble(Cols.size()*2);
    for(size_t c=0;c<Cols.size();c++)
        for(std::complex<double> v: Cols[c]){
            cDouble[2*c  ].push_back(v.real());
            cDouble[2*c+1].push_back(v.imag());
        }
    writeCols(cDouble,Digits);
}

void AsciiFile::writeCols(const std::vector<std::vector<double>> &Cols, unsigned int Digits) const {
    if(Cols.size()==0)return;

    ofstream stream;
    stream.open(_name.c_str(),std::ios_base::out|std::ios_base::app);
    if(not stream.is_open())ABORT("could not open input file '"+_name+"'");
    size_t maxRow=Cols[0].size();
    for(size_t l=1;l<Cols.size();l++){
        maxRow=std::max(maxRow,Cols[l].size());
        if(_blank<Cols.size() and Cols[l-1].size()!=Cols[l].size())
            ABORT("2d data storage only for all columns equal length");
    }

    for(size_t k=0;k<maxRow;k++){
        if(_blank<Cols.size() and k>0 and Cols[_blank][k-1]!=Cols[_blank][k])stream<<" "<<endl;
        for(unsigned int l=0;l<Cols.size();l++){
            if(Cols[l].size()>k)stream<<setprecision(Digits)<<Cols[l][k];
            if(l+1<Cols.size())stream<<_sep;
        }
        stream<<endl;
    }
}
void AsciiFile::copy(string NewFile, bool Overwrite){
    if(not Overwrite and folder::exists(NewFile))ABORT("for now, always use with Overwrite=true");
    ofstream out(NewFile.c_str());
    copy(out);
    out.close();
}
void AsciiFile::copy(ostream &Stream){
    ifstream inp(name().c_str());
    string line;
    while(getline(inp,line))Stream<<line<<endl;
}
void AsciiFile::writeRow(const std::vector<double> & Row, unsigned int Digits) const {
    ofstream stream;
    stream.open(_name.c_str(),std::ios_base::out|std::ios_base::app);
    if(not stream.is_open())ABORT("could not open input file '"+_name+"'");
    stream<<setprecision(Digits);
    for(size_t k=0;k<Row.size();k++){
        stream<<setw(Digits+7)<<Row[k];
        if(k<Row.size()-1)stream<<setw(_sep.length())<<_sep;
    }
    stream<<endl<<flush;
}
void AsciiFile::writeBlankRow(int Count) const {
    ofstream stream;
    stream.open(_name.c_str(),std::ios_base::out|std::ios_base::app);
    if(not stream.is_open())ABORT("could not open input file '"+_name+"'");
    for(int k=0;k<Count;k++)stream<<endl;
}

bool AsciiFile::empty() const{
    ofstream stream;
    stream.open(_name.c_str(),std::ios_base::out|std::ios_base::app);
    return stream.tellp()==0;
}

void AsciiFile::readCols(std::vector<Eigen::MatrixXd> &Cols, std::vector<unsigned int> NCols, const std::string Sep, bool SetsRowWise){
    ifstream stream(_name.c_str(),std::ios::in);
    if(not stream.is_open())ABORT("could not open input file '"+_name+"'");

    vector<int>siz;
    vector<vector<double> > cols;
    readCols(cols,NCols,Sep,siz);
    for(size_t k=0;k<cols.size();k++){
        if(SetsRowWise)Cols.push_back(Eigen::Map<Eigen::MatrixXd>(cols[k].data(),siz[0],cols[k].size()/siz[0]).transpose());
        else           Cols.push_back(Eigen::Map<Eigen::MatrixXd>(cols[k].data(),siz[0],cols[k].size()/siz[0]));
    }
}

void AsciiFile::readCols(std::vector<std::vector<double> > &Cols, std::vector<unsigned int> NCols, const string Sep, int DataSet){
    vector<int>siz,set(1,DataSet);
    if(DataSet==-1)set.clear();
    readCols(Cols,NCols,Sep,siz,set);
}

bool AsciiFile::addCols(const std::vector<std::vector<double> > & Cols, std::vector<size_t> Match){
    std::vector<std::vector<double>> cols;
    readCols(cols);
    if(cols.front().size()!=Cols.front().size())return false;
    if(Match.size()==2){
        if(cols[Match[0]].size()!=cols[Match[1]].size())return false;
        double vMax(DBL_MAX),vMin(-DBL_MAX);
        for(auto v: cols[Match[0]]){
            vMin=std::min(vMin,v);
            vMax=std::max(vMax,v);
        }
        for(size_t k=0;k<cols[Match[0]].size();k++){
            if(abs(cols[Match[0]][k]-Cols[Match[1]][k])>1.e-12*(vMax-vMin))return false;
        }
    }
    std::vector<std::string>comm;
    readComments(comm);
    writeComments(comm);
    for(size_t k=0;k<Cols.size();k++)
        if((Match.size()==2 and k!=Match[1]) or Match.size()==0)cols.push_back(Cols[k]);
    writeCols(cols);
    return true;
}

void AsciiFile::readCols(std::vector<std::vector<double> > &Cols, std::vector<unsigned int> NCols, const string Sep,
                         std::vector<int> &DataSize, const std::vector<int> DataSet)
{
    if(Sep.length()!=1)DEVABORT("colum separator must be length=1, is ->"+Sep+"<-");
    ifstream stream(_name.c_str(),std::ios::in);
    if(not stream.is_open())ABORT("could not open input file '"+_name+"'");

    string line;
    vector<string> item,seps;
    unsigned int nIn;

    int set=0;
    DataSize.assign(1,0);
    while(getline(stream,line))
    {
        if(line.find_first_of(comment)==0)continue; // skip comment lines

        if(line==""){
            set++;
            if(DataSet.size()!=0 and find(DataSet.begin(),DataSet.end(),set)==DataSet.end())continue;
            DataSize.push_back(0);
        }

        else{
            if(DataSet.size()==0 or find(DataSet.begin(),DataSet.end(),set)!=DataSet.end()){
                tools::splitString(tools::lcropString(line),Sep,item,seps);
                nIn=item.size();
                if(Cols.size()==0){
                    if(NCols.size()==0)
                        for(unsigned int k=0;k<item.size();k++)NCols.push_back(k);
                    if(item.size()==0 or tools::anyElement(NCols,tools::larger,(unsigned int)item.size()-1))
                        ABORT("too few items in line, NCols="+tools::str(NCols,3," ,")+"\n"+line);
                    if(NCols.size()==0)return;
                    Cols.resize(NCols.size());
                }

                if(item.size()!=nIn)ABORT("number of items in row changed");
                DataSize.back()++;
                for(unsigned int k=0;k<NCols.size();k++){
                    Cols[k].push_back(tools::string_to_double(item[NCols[k]]));
                }
            }
        }
    }
}

void AsciiFile::readCol(std::vector<double> &Col, unsigned int Pos){
    vector<vector<double> > cols;
    vector<unsigned int> pos(1,Pos);
    readCols(cols,pos);
    Col=cols[0];
}
