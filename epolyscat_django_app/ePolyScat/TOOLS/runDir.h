#ifndef RUNDIR_H
#define RUNDIR_H

#include <string>

/// manages tRecX's run-directory names
///
/// typical format NameOfInput/0123/
class RunDir
{
    std::string _root; // run root, includs trailing /
    int _digits;
    int _cur;
    int digits(); ///< determine numer of digits as actually used, _digits if no created yet
public:
    RunDir():_digits(4),_cur(0){}
    RunDir(std::string Root, int Digits=0):RunDir(){updateRoot(Root);updateDigits(Digits);}

    void updateRoot(std::string Root); ///< reset Root name, if changes, determine number of digits
    void updateDigits(int Digits); ///< reset number of digits for run directory names
    int cur() const {return _cur;} ///< current run number

    /// current directory name, format Root/Count/ (includes trailing /)
    /// ="" count at pow(10,_digits)
    std::string dir() const;
    RunDir& nextUnused(); ///< advance to next unused Root/Count directory
};


#endif // RUNDIR_H
