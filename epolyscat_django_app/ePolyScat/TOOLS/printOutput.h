// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef PRINTOUTPUT_H
#define PRINTOUTPUT_H
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "timer.h"
#include "str.h"
#include "mpiWrapper.h"

class MultiStream;

/// \ingroup IO
/// \brief unified formated output for project
///
/// output is by default to screen, can be duplicated to file, and can be suppressed
/// by specifying the "-noScreen" flag
class PrintOutput
{
public:

    static const std::string outExtension;   ///< default extension to use for output files
    static void set(std::string File);       ///< direct output to File (for direct "cout" to cout)
    static bool noScreenFlag();                 ///< check for -noScreen flag
    static void off(std::string Which="both");            ///< switch off all output from PrintOutput
    static void header(std::string File=""); ///< write header to File, File=="": do not change current
    static void title(std::string Title);    ///< write a title over group of outputs
    static void subTitle(std::string Sub);   ///< write subtitle over group of outputs
    static void message(std::string Mess, std::ostream *Ostr=0, bool Final=false, int Limit=1); ///< send a message to output
    static void matrix(const Eigen::MatrixXd & Mat,int Digits=3); ///< table-like print of matrix, entries<10^(-Digits) suppressed
    static void terminationMessages(std::ostream *Ostr=0); ///< reapeat messages listed as "Final=true"
    static void verbatim(std::string Value,std::ostream *Ostr=0); ///< put into output stream, formatting completely in Value
    static void progressStart(std::string Mess); ///< append Mess to line, flush
    static void progressStop(std::string Mess="",std::ostream *Ostr=0); ///< append Mess to line, flush
    static void progress(std::string Mess,std::ostream *Ostr=0); ///< append Mess to line, flush
    /// send a message to output (Limit<0: warn only after |Limit| calls, Explanations will be printed into file)
    static void warning(std::string Mess, int Limit=50, std::ostream *Ostr=0, const std::string Explanations="");
    static void warningList(); ///< send a message to output
    static void paragraph();                                    ///< end group of outputs with empty line
    static void newLine();                                      ///< end group, start new line
    static void lineItem(std::string Name,double      Value,std::string Subtitle="",unsigned int Digits=4); ///< add double item to line of output
    static void lineItem(std::string Name,int         Value,std::string Subtitle=""); ///< add int item to line of output
    static void lineItem(std::string Name,std::string Value,std::string Subtitle=""); ///< add string item to line of output
    static void tableColumn(std::string Name,std::vector<double>      Values,std::string Subtitle="");
    static void tableColumn(std::string Name,std::vector<int>         Values,std::string Subtitle="");
    static void tableColumn(std::string Name,std::vector<std::string> Values,std::string Subtitle="");
    static void tableColumn(std::string Name,std::string Subtitle=""); ///< start new column
    static void newRow();                   ///<open a new row in table
    static void rowItem(double Value, unsigned int Width=4); ///< append double item to current row
    static void rowItem(       int Value);  ///< append int item to current row
    static void rowItem(unsigned int Value);///< append unsigned int item to current row
    static void rowItem(std::string Value); ///< append string item to current row
    static void flush(); ///< write current output, without terminating
    static void updateTable(); ///< append new outputs to table
    static void end(); ///< flush and explicitly terminate any set of outputs (line,table...)
    static std::map<std::string,unsigned int>warningCount;
    static std::map<std::string,unsigned int>messageCount;
    static std::vector<std::string>finalMessage;

    static std::ostream * output(){if(sOut!=0)return sOut; return &std::cout;} ///< return output stream (file or cout)

    static void outputLevel(std::string Level);
    static void timerWrite(); ///< write timer out to output stream(s)
    static std::string level;

    static void DEVmessage(std::string Mess, std::ostream *Ostr=0, bool Final=false);///< message for developers
    static void DEVwarning(std::string Mess, int Limit=50, std::ostream *Ostr=0);///< warning for developers

    static std::string indent(){return _indent;}
    static std::string columnSeparator(){return _columnSeparator;}
private:
    static MultiStream * out;
    static std::ofstream * sOut;
    static std::ofstream * sMon;
    static std::string fileOut;
    static std::string monMessage;

    static unsigned int lengthLine;
    static unsigned int namesTab;
    static std::string _indent;
    static std::string titleSeparator;
    static std::string itemSeparator;
    static std::string nameSeparator;
    static std::string _columnSeparator;
    static std::string namesEnd;
    static std::string namesLine,valuesLine;

    static std::string _off;

};

#endif // PRINTOUTPUT_H
