#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H
#include <stdlib.h>
#include <cfloat>
#include <cmath>
/// solver for x: Op x = b
///
/// with A: V->V an inertible linear map
template<class Oper, class V>
class LinearSolver{
protected:
    const Oper* _op;
    static constexpr double defaultEps=DBL_MIN*100;
    double _eps;
    mutable int _nSolves;
public:
    LinearSolver(const Oper* Op):_op(Op),_eps(defaultEps),_nSolves(0){}

    /// possible setup step (e.g. after update of Op)
    virtual void compute(const Oper *Op)=0;

    /// set acceptable L2-error
    LinearSolver& withTolerance(double Eps){_eps=Eps; return *this;}

    /// return acceptable L2-error
    double tolerance() const {return _eps;}

    /// number of calls to solve()
    int nSolves() const {return _nSolves;}

    /// returns x= Op^-1 B
    virtual V& solve(const V& B) const =0;

    bool test(V& Vtest /** set to non-zero Vtest */){
        double nrmsq=std::abs(Vtest.dotProduct(Vtest));
        if(nrmsq<tolerance()*tolerance())std::abort();
        Vtest*=1./std::sqrt(nrmsq);
        V res(Vtest);
        _op->apply(1.,Vtest,0.,res);
        Vtest-=solve(res);
        return Vtest.dotProduct(Vtest)<tolerance()*tolerance();
    }
};


#endif // LINEARSOLVER_H
