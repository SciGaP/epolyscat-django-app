#ifndef PROBLEMRESULT_H
#define PROBLEMRESULT_H

#include "problemNode.h"

class ProblemParameters;

/// info about a specific parameter point
///
///<br> errors at a set of specific parameter values
///<br> estimates in for parameters in the vicinity
///<br> approximate behavior of the error near that point
///<br> real times and total CPU times at machine and MPI nodes
class ProblemResult: public ProblemNode
{
protected:
    double _cpu;
    std::string _host;
    int _nproc;
    std::string _dir;

    struct res{
        std::vector<std::string> _def; // parameter. In tRecX: category:name[line]
        int _ref; // _ref... reference parameter value (must match one of the _pars)
        std::string _err; // _err...algebra string for error estimate near _ref
        std::vector<double> _pars,_errs; // all _pars and matching error estimates of reference calculation
        res(std::string Line);
        std::string parnam() const {return _def[0]+":"+_def[1]+(_def[2]=="0"?"":"["+_def[2]+"]");}
        std::vector<std::string>parsErrs() const {
            std::vector<std::string> res;
            for(size_t k=0;k<_pars.size();k++)
                res.push_back(tools::joinString({tools::str(_pars[k],3),tools::str(_errs[k],3)},":"));
            return res;
        }
        std::string line() const;
        int ref() const {return std::find(_pars.begin(),_pars.end(),_ref)-_pars.begin();}
    };
    std::vector<res>_res; // list of all parameters in varied for error estimate

    double errTotal() const;
public:
    bool hasConvergence() const; // true if convergence information is available
    ProblemResult():ProblemResult("Result"){}
    std::string definition() const {
        std::vector<std::string> parts({tools::str(_cpu,3),_host,tools::str(_nproc),_dir});
        for(auto r: _res)parts.push_back(r.line());
        return definitionJoin(parts);
    }
    ProblemResult(std::string Definition)
        :ProblemNode("Result"){
        auto parts=definitionSplit(Definition);
        if(parts.size()==0)return;

        _cpu=tools::string_to_double(parts[0]);
        _host=parts[1];
        _nproc=tools::string_to_int(parts[2]);
        _dir=parts[3];
        for(size_t k=4;k<parts.size();k++)_res.push_back(res(parts[k]));
   }

    void addToInput(std::ofstream &Stream) const{Stream<<"";};
    void addToHeader(std::ofstream &Stream) const;
    std::string str() const{return "Accuracy: "+tools::str(errTotal(),3)+" at CPU "+tools::str(_cpu*_nproc);}
};

#endif // PROBLEMRESULT_H
