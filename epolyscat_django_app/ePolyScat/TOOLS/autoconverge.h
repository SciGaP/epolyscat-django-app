#ifndef AUTOCONVERGE_H
#define AUTOCONVERGE_H


class ReadInput;
class Farm;
class Metric;

#include "tree.h"
#include "convergenceNode.h"
#include "inputGenerator.h"

#include "autoTarget.h"
class TargetError;

/// automatically control convergence
class Autoconverge: public Tree<Autoconverge>, public InputGenerator
{
    static std::string estimate;

    //----- variables that should only be at root ----------------
    std::string _rootInp; // root of all input files
    std::string _reference; // reference input file at present iteration

    double _accuracyTarget; // desired accuracy of target measure
    std::shared_ptr<AutoTarget> _target; // definition of target data for convergence
    std::shared_ptr<Metric> _metric; // how to compute distance between target columns
    //------- access to the above variables from below --------------------------

    std::shared_ptr<ConvergenceNode> _conNode;

    int parameterRow() const; // row where _cat and _name appear

    /// true if converged wrt present parameter
    bool converged(){return false;}///<dummy for now
    bool addNextParameter(size_t Priority); // extend _pars

    /// compute the distance between target data in DirA and DirB
    double distance(std::string DirA, std::string DirB) const;

    /// find completed runs and compute delta
    void getDeltas(std::vector<double> &CPars, std::vector<double> &Deltas, std::vector<double> &Gammas, std::vector<double> &Alphas) const;

    std::string writeConvergence(); ///< write convergence data to file, return name of file

    std:: string strNode(int Level=Tree_defaultKind) const;
    void writeResults(std::ofstream &ResultStream) const;
    const ConvergenceNode * node() const {return _conNode.get();}
public:
    static const std::string macro;
    static std::shared_ptr<Autoconverge> read(ReadInput& Inp){return std::shared_ptr<Autoconverge>(new Autoconverge(Inp));}
    Autoconverge(ReadInput & Inp, std::vector<Autoconverge *> Path={});

    std::string nextPriority(size_t Priority); ///< return next inputs at priority levels > Priority, 10 - highest

    void print(); ///< print status
    bool runCompleted(std::string RunDir) const;
    std::string nextInput();
    int nInps() const;
    /// actually computed parameters and error estimates
    void parametersAndEstimates(std::vector<double> &Pars, std::vector<double> &Estim) const;
    /// write convergence resutls summary to file, return file name
    std::string writeResults() const;

    double accuracyTarget() const {return root()->_accuracyTarget;}
    std::shared_ptr<AutoTarget> target() const {return root()->_target;}
    std::shared_ptr<Metric> metric() const {return root()->_metric;}
    std::string reference() const {return root()->_reference;}
    std::string referenceRoot() const  {return root()->_rootInp;}

};

#endif // AUTOCONVERGE_H
