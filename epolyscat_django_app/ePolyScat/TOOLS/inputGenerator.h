#ifndef INPUTGENERATOR_H
#define INPUTGENERATOR_H

#include <memory>
#include <vector>
class ReadInput;

#include "reuseRun.h"
#include "runDir.h"

class InputGenerator{
protected:
    // control re-use of previous runs
    std::shared_ptr<ReuseRun> _reuse;
    RunDir _runDir;

    std::string writeInpFile(std::string Root, const std::vector<std::string> &InputLines);
    std::vector<std::string> readInpFile(std::string InpFile);
public:
    virtual void print(){} ///< print generator status
    virtual bool runCompleted(std::string RunDir) const {RunDir=""; return false;} ///< overwrite with completion criteria in derived classes
    static std::shared_ptr<InputGenerator> factory(int argc, char* argv[]);
    virtual std::string nextInput() =0; ///< return path to next input
    virtual int nInps() const =0; ///< expected number of concurrent inputs (for partioning farm into flocks)
};

class InputSingle: public InputGenerator{
    mutable std::string _nxt;
public:
    /// present a single input
    InputSingle(){};
    InputSingle(int argc, char* argv[]);
    std::string nextInput(){std::string res(_nxt);_nxt="";return res;}
    int nInps() const {return int(_nxt!="");}
};

class InputRanges: public InputGenerator{
    std::string _root;
    std::vector<std::vector<std::string>> _inpLines;
    mutable int _current;
public:
    /// generates inputs by looping through input item values
    ///
    /// loops can be specified in one of the following forms
    /// <br> RANGE{a,b,c,d} for all elements of the list
    /// <br> RANGE{a:b:n}   for n equidistant floats from first a to last b
    InputRanges(int argc, char* argv[]);
    std::string nextInput();
    int nInps() const;
};


#endif // INPUTGENERATOR_H
