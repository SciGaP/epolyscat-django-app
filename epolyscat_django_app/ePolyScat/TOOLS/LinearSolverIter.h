#ifndef LINEARSOLVERITER_H
#define LINEARSOLVERITER_H

#include "LinearSolver.h"
#include <vector>

template<class Oper, class V>
class LinearSolverIter:public LinearSolver<Oper,V>{
protected:
    const int _maxIter;
    mutable int _nIter;
    mutable int _totIter;
    mutable int _nApplys;
    const V* _guess;
    mutable std::vector<double> _residues;
public:
    LinearSolverIter(const Oper* Op,int MaxIter):LinearSolver<Oper,V>(Op),
        _maxIter(MaxIter),_totIter(0),_guess(0){}

    /// set a starting vector for the iteration
    LinearSolverIter& withInitialVector(const V* Initial){_guess=Initial; return *this;}

    /// prepare as needed
    void compute(const Oper *Op)=0;
    virtual V& solve(const V&) const=0;

    /// number of iterations
    int nIter() const {return _nIter;}
    int totIter() const {return _totIter;}
    int nApplys() const {return _nApplys;}
    bool converged() const {return _nIter<_maxIter;}
    const std::vector<double> & residues() const {return _residues;}

};


#endif // LINEARSOLVERITER_H
