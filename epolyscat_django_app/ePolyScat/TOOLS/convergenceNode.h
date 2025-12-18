#ifndef CONVERGENCENODE_H
#define CONVERGENCENODE_H


class ReadInput;
class Farm;
class Metric;

#include "tree.h"
#include "autoTarget.h"
#include "targetError.h"

class TargetError;
class Autoconverge;

/// convergence wrt to single parameter: estimate errors, propose next parameter
class ConvergenceNode
{
    static std::string estimate;
    // note: we organize this as a tree for possibly more complicated hierarchies
    // each node handles one parameter

    const std::shared_ptr<AutoTarget> target() const;
    std::shared_ptr<Metric> metric() const;
    std::string reference() const;
    std::string referenceRoot() const;
    size_t _size;


protected:
    const Autoconverge* _root;

    std::string _cat,_nam;
    double _parStep,_parMin,_parMax; // initial step size and limits for parameter
    mutable std::vector<double> _pars; // parameter values at present node
    mutable std::vector<std::string> _inps; // input files where present parameter varies
    int _row,_col; // row and column of present parameter in input matrix
    mutable std::shared_ptr<TargetError> _targetError;
public:

    std::string cat() const {return _cat;}
    std::string nam() const {return _nam;}
    double parStep() const {return _parStep;}
    double parMin() const {return _parMin;}
    double parMax() const {return _parMax;} // initial step size and limits for parameter
    std::vector<double> & pars() const {return _pars;} // parameter values at present node
    std::vector<std::string> & inps() const {return _inps;} // input files where present parameter varies
    int row() const {return _row;}
    int col() const {return _col;} // row and column of present parameter in input matrix
    std::shared_ptr<TargetError> targetError() const {return _targetError;}
    size_t size() const {return _size;}

protected:
    int parameterRow() const; // row where _cat and _name appear

    /// true if converged wrt present parameter
    bool converged(){return false;}///<dummy for now

    /// compute the distance between target data in DirA and DirB
    double distance(std::string DirA, std::string DirB) const;


    std::string writeConvergence(); ///< write convergence data to file, return name of file


    /// input as a matrix of strings, return in _row, _col location of parameter in matrix
    std::vector<std::vector<std::string>> inputRowCol(ReadInput& Inp, std::string AtLine, int &Row, int &Col);
    void writeResults(std::ofstream &ResultStream) const;


    mutable std::vector<std::vector<std::string>> _inputRows; // input as a list of input rows
public:
    /// most recent error fit
    std::shared_ptr<TargetError> & errorFunction() {return _targetError;}
    std::shared_ptr<TargetError> & errorFunction() const {return _targetError;}


    /// find completed runs and compute delta
    void getDeltas(std::vector<double> &CPars, std::vector<double> &Deltas) const;
    /// get estimates of errors (depends on errorFunction)
    void parametersAndErrors(std::vector<double> &Pars, std::vector<double> &Estims) const;

    static const std::string macro;
    ConvergenceNode():_size(0),_root(0),_row(-1),_col(-1){}
    ConvergenceNode(ReadInput & Inp, const Autoconverge* Root, int line);

    /// algorithm to suggest next parameter for convergence
    virtual double nextParameter(std::vector<double> &Pars, std::vector<double> &EstimErrs, double &ParAtTarget) const;

    std::vector<std::vector<std::string>> & inputRows() const {return _inputRows;}

    /// write convergence resutls summary to file, return file name
    std::string writeResults() const;

    std::string strNode(int Level) const;

};

#endif // CONVERGENCENODE_H
