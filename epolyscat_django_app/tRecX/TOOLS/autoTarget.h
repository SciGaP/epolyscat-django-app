#ifndef AUTOTARGET_H
#define AUTOTARGET_H

#include "asciiFile.h"
#include "convergenceTarget.h"

class AutoTarget:public ConvergenceTarget{
    std::string _fileKind; // file type
    mutable std::string _rows;  // code for rows to be selected
    std::vector<unsigned int> _cols; // column numbers on file
    std::vector<std::string> _defaults; // defaults for various file types

    std::vector<double> rowPars() const;
    bool inTarget(int K, std::vector<double> RefCol, std::vector<double> Pars);
    std::string file() const {return _fileKind;};

public:
    ///@brief Target observable for tRecX computation (column-wise file)
    /// Def=file[rows][col1,col2,...]
    ///<br> rows=all or rowfront:rowback or minVal,maxVal,xAxis   where xAxis defaults to 0
    ///<br>  rows or both specifications can be omitted
    AutoTarget(std::string Def, std::vector<std::string> Defaults):ConvergenceTarget("AutoTarget"){_defaults=Defaults;_def=Def;}
    std::vector<std::vector<double>> cols(std::string Dir); ///< columns from Dir/file()
    virtual std::string str() const;
    std::string definition() const {return name()+": "+str();}

};

#endif // AUTOTARGET_H
