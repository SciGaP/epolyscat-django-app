// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef READINPUT_H
#define READINPUT_H
#include <vector>
#include <set>
#include <string>
#include "folder.h" // for handling folders
#include <climits>
#include <stdio.h>
#include <map>
#include "inputitem.h"

#ifdef _WIN32
// for alternative tokens for windows ( e.g. and, or, not, ... ) 
#include <iso646.h>
#endif


/** @defgroup IO Input and output
 *  @ingroup Tools
  * \brief read control input, standard structured output, open, write and read standard data files, etc.
  * @{
  */
/// @brief Universal input
///
/// all human-readable input should be covered by this class
///
/// see static Test for usage examples
///
/// input data after a marker line of the form
/// "Category: name1, name2, name3"
/// at the corresponding column and at specified line line (default=first line)
/// all values can be overwritten by command line flags
/// default flag -Category:name1=value, flag can be specified explicitly
class ReadInput{
    friend class ReadInputRange; // needs access to all items
    std::vector<std::string> commandLine;
    // given an input file, construct output directory and root names
    void rootAndDir(std::string File);
protected:
    size_t posExactName(std::string Line, std::string Name); // find postion of exact name (item between commas)
    virtual std::string readValue(const std::string Category, const std::string Name, const std::string Default, const std::string Docu,
                                  unsigned int Line, std::string Flag, std::string Allow);
    bool inCommandLine(std::string Flag,std::string DefaultFlag,unsigned int Line,std::string & valString);
    std::string strAdmissibleInputs(std::string Docu) const; ///< line-break separated list of admissible inputs (form Item::docu)
    void error(std::string Message, unsigned int Item, unsigned int Item2=INT_MAX) const;
    void flagError(std::string Flag) const;
    void parseInput(std::string Line,std::string & Category, std::vector<std::string > & Name);

    void inputMacros(); ///< substitute input macros

    // all data can be pointed to
    std::vector<InputItem> inputTable;  ///< keeps track of all input items in present ReadInput
    std::string unitSystem; ///< if input units but no output units are given, will be converted to this unitSystem
    std::string mainProgram; ///< name of the main program
    std::string _root;
    mutable std::string _outDir;
    std::string _outPrefix;
    std::string subDir; // subDir replaces _outPrefix
    std::vector<std::string> com; //<command line input strings

    std::vector<std::string> inpLines;
    std::vector<std::string> allLines;

    bool _writeLinp;

    InputItem & findMatching(InputItem &Cur);
    bool isDuplicate(InputItem &Cur){return &Cur!=&findMatching(Cur);}
    bool isInfinity(std::string Inf){return Inf.find("Infty")<2 or Inf.find("inf")<2;}
    bool   string_to_bool(  const std::string &);
    int    string_to_int(   const std::string &);
    double string_to_double(const std::string &);

    void construct(const std::string & File, int NCom, char * Comp[], bool FromCommand);
    static InputItem dummyItem;
    InputItem & readItem(const std::string Category, const std::string Name, const std::string Default, const std::string Docu,
                         unsigned int Line, std::string Flag, std::string Allow);

    static unsigned int mpiMaster;

    // return Value if in allowed list
    double      allowedNumber(std::string Value,InputItem& It);
    std::string allowedString(std::string Value,InputItem& It);
public:

    static const std::string flagFound;  ///< indicate that flag was found, but no value specified
    static const std::string notFound;   ///< indicate that Category:name was not specified
    static const std::string empty;      ///< indicate that Category:name is in input file, but entry is blank

    // string characters
    static const std::string itemSeparator;     ///< items are separated by these characters
    static const std::string categoryTerminator;///< categories are terminated by this character
    static const std::string whiteSpace; ///< list of valid white space characters
    static const std::string comments;  ///< list of valid comment characters
    static const std::string quote;     ///< admissible quote characters
    static const std::string leftBracket;  ///< all known left brackets (matches pairwise with rightBracket)
    static const std::string rightBracket; ///< all known left brackets (matches pairwise with leftBracket)
    static const std::string forbidden; ///< forbidden input characters

    // file extensions:
    static const std::string inputExtension; ///< input files must end in this extension
    static const std::string inputCopy; ///< name of input copy in data directory
    static const std::string inputList; ///< file extension for list of all inputs (mostly for machine reading)
    static const std::string docExtension; ///< extension for input documentation file

    // special input item strings
    static const std::string flagOnly;   ///< indicates that -flag was specified without value
    static const std::string noDefault;  ///< cannot leave unspecified
    static const std::string noItem;              ///< do not insert into item list
    static const std::string doNotAddToItemTable; ///< do not insert into item list
    static const std::string anyName; ///< any name of category
    static const std::string macro;   ///< string to indicate a macro

    std::map<std::string,std::string> readList; ///< keeps track of multiple reads
    static void mpiSetMaster(unsigned int Master){mpiMaster=Master;}

    static std::string getInputFile(const std::string & File /**< input file name */,
                                    int NCom=0, /** 0...no flag input, -1...no flag and no warning*/
                                    char * Com[]=0, bool FromCommand=true);
    static std::string nextRunDir(std::string Root); ///< create new Root/0123 type directory and return name

    ~ReadInput(){}
    ReadInput():_writeLinp(true){}
    ReadInput(const std::string & File /**< input file name */,
              int NCom=0, /** 0...no flag input, -1...no flag and no warning*/
              char * Com[]=0, bool FromCommand=true); ///< create a new (local) input
    ReadInput(int NCom, char * Com[]); ///< input file = first command line argument

    ReadInput & suppressLinp(){_writeLinp=false; return *this;}

    static ReadInput main;              ///< main input (globally available)
    static void openMain(const std::string & File /**< input file name */,
                         int Ncom=0, char *Com[]=0, bool FromCommand=true){ReadInput::main=ReadInput(File,Ncom,Com,FromCommand);} ///< open the "main" input
    static void Test(int=0, char *Com[]=0); ///< usage examples

    void outPrefix(std::string Prefix){_outPrefix=Prefix+"/";} ///< Prefix for output()
    std::string output(){return outputTopDir()+subDir;} ///< return outputDirectory+Prefix+"_"
    std::string outputTopDir() const;/// generate an output directory return name, derive name from input file
    std::string outputFile();//{return output()+PrintOutput::outExtension;}
    void outSubdir(std::string Name);

    /// close input file, check for illegal items, create output directory, copy input (flags -DEBUG.... will not be checked)
    void finish();

    /// display ReadInput
    void show() const;
    static void showHelp(std::vector<std::string> CommandLine, std::string Category="");
    virtual void writeDoc(std::ostream *Doc=0, std::vector<std::string> AddLines={}); ///< write input item documentation to file
    void writeLinp() const;
    virtual std::string root(){return _root;} ///< root to all runs in with same input name
    virtual std::string file(){return _root+inputExtension;} ///< input file name
    virtual std::string docFile(){return mainProgram+docExtension;} ///< input file name
    void texdocuCategoryAdd(std::string Category, std::string Sort, std::string TexString, std::string Tutorials="");

    /// stop and send message if command has become obsolete
    void obsolete(std::string Category,std::string Name, std::string Message,unsigned int Line=1,std::string Flag="");

    /// @brief retrieve string value of input item
    ///
    /// input parameters: see version for "double"<br>
    /// leading and trailing blanks will be cropped<br>
    /// outermost quotes will be stripped, if quotes are required, include within pair of other type: " or '
    InputItem & read(std::string Category, std::string Name, std::string & Value, std::string Default, std::string Docu, unsigned int Line=1, std::string Flag="", std::string Allow="");

    /// @brief read double values, automatic unit conversion
    ///
    /// Example:
    ///
    /// inp.read("Field","Intensity(W/cm2)",inten,"0.01","some intensity",2) <br>
    /// inp.read("Field","wavelength(nm)",lambda,"100.","some wave length") <br>
    /// inp.read("Field","frequency",freq,"1.","some frequecy") <br>
    /// inp.read("Field","shape",shape,"1",","can have values: sin2...sine-square, gauss...gaussian") <br>
    /// <br>Field: shape, wavelength(nm), Intensity(W/cm2), frequency(au) <br>
    ///    gauss,     800   ,            0 ,   2/pi au <br>
    ///    sin2,      400   ,         1 au <br>
    /// <br>
    /// results in:<br>
    /// lambda=800: no units given, assumed to be in the specified units <br>
    /// inten =3.509e16: 2nd item in 2nd line, atomic units "au" converted to "W/cm2" <br>
    /// freq =1/3.1415...(exact): one can use algebraic expresion and certain predifined constants <br>
    /// shape = gauss: the values following : and terminated by ... are the only admissible inputs, escape by starting docu with {allow all}
    InputItem & read(std::string Category, std::string Name, ///< free-formated string, can specify return units, e.g. wavelength(nm), Intens(W/cm2),etc.
                     double & Value, ///< return value (optionally converted to units as specified in Name
                     std::string Default, ///< default input (string!)
                     std::string Docu,    ///< brief description of input item
                     unsigned int Line=1, ///< line after Name-line from where to read value
                     std::string Flag="", ///< command line flag (default is Category:Name)
                     std::string Allow="" ///< comma-separated list of allowed input values
            );
    /// input parameters: see version for "double"
    InputItem & read(std::string Category, std::string Name,        bool & Value, std::string Default, std::string Docu, unsigned int Line=1, std::string Flag="");
    /// input parameters: see version for "double"
    InputItem & read(std::string Category,std::string Name,         int & Value,std::string Default,std::string Docu,unsigned int Line=1,std::string Flag="",std::string Allow="");
    /// input parameters: see version for "double"
    InputItem & read(std::string Category,std::string Name,unsigned int & Value,std::string Default,std::string Docu,unsigned int Line=1,std::string Flag="",std::string Allow="");
    /// read a comma-separated list into vector
    /// input parameters: see version for "double"
    InputItem & read(std::string Category,std::string Name,std::vector<std::string> & Value,std::string Default,std::string Docu,unsigned int Line=1,std::string Flag="",std::string Allow="");
    /// input parameters: see version for "double"
    InputItem & read(std::string Category,std::string Name,std::vector<int>         & Value,std::string Default,std::string Docu,unsigned int Line=1,std::string Flag="",std::string Allow="");
    /// input parameters: see version for "double"
    InputItem & read(std::string Category,std::string Name,std::vector<double>      & Value,std::string Default,std::string Docu,unsigned int Line=1,std::string Flag="",std::string Allow="");

    double readDouble(std::string Category,std::string Name,std::string Default,std::string Docu,unsigned int Line=1,std::string Flag="",std::string Allow="");
    int    readInt   (std::string Category,std::string Name,std::string Default,std::string Docu,unsigned int Line=1,std::string Flag="",std::string Allow="");
    bool   readBool  (std::string Category,std::string Name,std::string Default,std::string Docu,unsigned int Line=1,std::string Flag="",std::string Allow="");


    /// flag only input (true if flag is specified)
    bool flag(std::string Name,std::string Docu){bool val;read("Flag",Name,val,flagOnly,Docu,0,Name);return val;}
    /// true if -Name[=value] is on commmand line
    bool hasCommandLineFlag(std::string Name);

    /// true if found in input file or command line
    bool found(std::string Category,std::string Name=anyName,std::string Flag=anyName);

    /// line number(s) of Category:name after Category, 0 for flag
    std::vector<size_t> foundAtLines(std::string Category, std::string Name);

    /// true if line is a syntactically valid Category line
    static bool isCategoryLine(std::string Line);

    /// list of all categories with format "Some_Name[specification]:"
    std::vector<std::string> categoryWithSpecification(std::string Category);

    /// list of all specs from entries Category:Name[spec]
    std::set<std::string> namesWithSpecification(std::string Category,std::string Name);

    /// true if Line'th line after first appearance of Category is empty
    bool endCategory(std::string Category,int Line);

    /// string at Line'th line after Category
    std::string lineAt(std::string Category, int Line=1);

    /// raw input lines (before any extension of macros
    std::vector<std::string> inputLines() const {return inpLines;}

    /// overwrite standard input by machine-generated list (TO BE IMPLEMENTED)
    void listMachine(std::string Category,std::string Name); ///< add item to list of machine-overwriteable items
    void readMachine(std::string,std::string Separator=","); ///< overwrite items by values in string (sequence must be as generated with listMachine)

    void setUnits(const std::string Units); ///< set default unitSystem for output
    void exclude(std::string Cat1, std::string Cat2, std::string Nam1=anyName, std::string Nam2=anyName); // define mutually exclusive input categories

    bool noInputFile() const; ///< true if no input file was given

    bool operator==(const ReadInput & Other) const;

};
/** @} */
#endif // READINPUT_H
