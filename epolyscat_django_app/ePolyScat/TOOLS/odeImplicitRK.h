// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef ODEIMPLICITRK_H
#define ODEIMPLICITRK_H

#include "odeStep.h"

#include <vector>
#include <complex>
#include <memory>
#include "abort.h"
#include "unistd.h"
#include "str.h"
#include "parameters.h"
#include "linSpaceMap.h"
#include "LinearSolver.h"
#include "LinearBiCGstab.h"
#include "MultiVector.h"
#include "butcherTableau.h"

#include "timer.h"
/// multi-vector operator for implicit solver
template<class Oper, class V>
class MultiOp: public LinSpaceMap<MultiVectorHilbert<V>> {
    Oper* _op;
    std::vector<std::vector<double>> _a;
    const std::vector<double> _c;
    double _time;
    double _step;
    const MultiVectorHilbert<V>* _model;
    std::vector<bool>_anyA;
    mutable size_t _opApply;
public:
    /// MultiOp(H,t,step) X: X(k) - step*H(t+step*c[k]){ Sum[l] a[k][l]*X(l) }
    MultiOp(Oper* OperOde,std::vector<std::vector<double>> A,std::vector<double> C,const MultiVectorHilbert<V>* Model)
        :_op(OperOde),_a(A),_c(C),_model(Model),_opApply(0){
        for(auto a: _a){
            _anyA.push_back(false);
            for(auto aa: a)_anyA.back()=_anyA.back() or aa!=0.;
        }
    }
    void update(double Time, const MultiVectorHilbert<V> * Dum){_time=Time;}
    void update(double Time, double Step){_time=Time,_step=Step;}

    /// Y(k)= A*{ X(k)- _step*H(t+_step*c[k])sum[l] a[k][l] X(l)) } + B*Y(k)
    void apply(std::complex<double> A, const MultiVectorHilbert<V> & X,std::complex<double> B, MultiVectorHilbert<V> & Y) const
    { // do only non-trivial applies
        if(A!=0.){
            for(size_t k=0;k<_a.size();k++){
                if(_anyA[k]){
                    _op->update(_time+_step*_c[k]);
                    _op->apply(-A*_step,X.reduceToSingle(_a[k]),B,Y(k));
                    _opApply++;
                }
                else
                    Y(k)*=B;
                Y(k).axpy(A,X(k));
            }
        }
        else Y*=B;
    }

    const MultiVectorHilbert<V> & lhsVector() const{return *_model;};
    const MultiVectorHilbert<V> & rhsVector() const{return *_model;};

    size_t opApplys() const {return _opApply;}
};


/// General Runge-Kutta time-stepper
///
/// default solver for implicit is BiCGstab
template<class Oper, class V>
class OdeImplicitRK: public OdeStep<Oper,V>{

    using OdeStep<Oper,V>::derOde; // C++? surprising it does not know where it is derived from

    std::vector<std::vector<double>> _butcherA;
    std::vector<double> _butcherB;
    std::vector<double> _butcherC;
    int _consistency;

    std::shared_ptr<MultiOp<Oper,V>> _multiOp;

    mutable std::shared_ptr<LinearSolver<MultiOp<Oper,V>, MultiVectorHilbert<V>> > _slv;
    LinearSolverIter<MultiOp<Oper,V>,MultiVectorHilbert<V>>*_iterSlv;
    double _solverEps;

    mutable MultiVectorHilbert<V> _vintern;
public:
    virtual ~OdeImplicitRK(){}

    /// SolverEps should be well below the desired step accuracy, but not excessivly small
    /// try e.g. SolverEps=accuracy/100
    OdeImplicitRK(Oper*Op,
                  double SolverEps /** tolerance */,
                  std::string Tableau)
        :OdeStep<Oper,V>("RK["+Tableau+"]",Op),_solverEps(SolverEps)
    {
        _butcherA=ButcherTableau(Tableau).a();
        _butcherB=ButcherTableau(Tableau).b();
        _butcherC=ButcherTableau(Tableau).c();
        _consistency=ButcherTableau(Tableau).consistency();

        _vintern=MultiVectorHilbert<V>(_butcherA.size(),derOde->lhsVector());
        _multiOp.reset(new MultiOp<Oper,V>(Op,_butcherA,_butcherC,&_vintern));
    };

    OdeImplicitRK<Oper,V> & withSolver(std::shared_ptr<MultiOp<Oper,V>> Slv){
        _slv=Slv;
        _iterSlv=dynamic_cast<LinearSolverIter<MultiOp<Oper,V>,MultiVectorHilbert<V>>*>(_slv.get());
        return *this;
    };

    OdeImplicitRK<Oper,V> & withSolverTolerance(double Eps){_solverEps=Eps;return *this;};


    /// single implicit RK step
    V &step(V &Vec, double Tstart, double  Tstep){

        if(not _slv){
            _slv.reset(new LinearBiCGstab<MultiOp<Oper,V>,MultiVectorHilbert<V>>(_multiOp.get(),10));
            _slv->withTolerance(_solverEps==0.?1.e-12:_solverEps);
            _iterSlv=dynamic_cast<LinearSolverIter<MultiOp<Oper,V>,MultiVectorHilbert<V>>*>(_slv.get());
        }

        // populate _vintern with H(t+c[k]*h) Vec
        for(size_t k=0;k<_vintern.mult();k++){
            derOde->update(Tstart+Tstep*_butcherC[k]);
            derOde->apply(1.,Vec,0.,_vintern(k));
        }

        // solve for K: [1-H(t,h)*A] K = H(t,h) Vec
        // solution is Vec->Vec+=sum[i] K(i)*b(i)
        _multiOp->update(Tstart,Tstep);
        _slv->compute(_multiOp.get()); // update solver as needed
        std::vector<double> hb(_butcherB);
        for(auto &b: hb)b*=Tstep;
        Vec+=_slv->solve(_vintern).reduceToSingle(hb);

        OdeStep<Oper,V>::nCallsStep=_multiOp->opApplys()-OdeStep<Oper,V>::nCalls;
        OdeStep<Oper,V>::nCalls=_multiOp->opApplys();
        return Vec;
    }

    /// return model vector argument of ODEstep
    const V & modelVec() const {return derOde->lhsVector();}

    unsigned int consistencyOrder() const {return _consistency;}
    double safetyFactor() const {return 0.1;}

    virtual std::string info() const {
        std::string s=OdeStep<Oper,V>::info();
        s+=(_iterSlv?", <iter>="+tools::str(double(_iterSlv->totIter())/_slv->nSolves(),3):"");
        return s;
    }

};

#endif // ODEIMPLICITRK_H
