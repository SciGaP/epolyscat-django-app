// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef ODECRANKNICOLSON_H
#define ODECRANKNICOLSON_H

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

#include "timer.h"


/// Crank-Nicolson time stepper y(2h)=(1-h Op[t+h])^-1 (1+h Op[t+h]) y(0)
///
/// default solver is BiCGstab
template<class Oper, class V>
class OdeCrankNicolson: public OdeStep<Oper,V>{
    using OdeStep<Oper,V>::derOde; // C++? surprising it does not know where it is derived from
    std::shared_ptr<V> z;
    /// implement (1-h Op)
    class One_hOp: public LinSpaceMap<V>{
        const Oper* _op;
        double _hHalf;
    public:
        One_hOp(const Oper* OperOde):_op(OperOde){}
        void update(double HHalf){_hHalf=HHalf;}
        void apply(std::complex<double> A, const V & X,std::complex<double> B, V & Y) const {
            _op->apply(-_hHalf*A,X,B,Y);
            Y.axpy(A,X);
        };
        const V & lhsVector() const{return _op->lhsVector();};
        const V & rhsVector() const{return _op->lhsVector();};

    };
    std::shared_ptr<One_hOp> _one_hOp;
    mutable std::shared_ptr<LinearSolver<Oper,V>> _slv;
    LinearSolverIter<Oper,V>*_iterSlv;
    double _solverEps;
public:
    virtual ~OdeCrankNicolson(){}

    /// SolverEps should be well below the desired step accuracy, but not excessivly small
    /// try e.g. SolverEps=accuracy/100
    OdeCrankNicolson(Oper*Op,double SolverEps /** tolerance */)
        :OdeStep<Oper,V>("CrankNicolson",Op),_solverEps(SolverEps)
    {
        _one_hOp.reset(new One_hOp(Op));
        z.reset(new V(derOde->lhsVector()));
    };

    OdeCrankNicolson<Oper,V> & withSolver(std::shared_ptr<LinearSolver<Oper,V>> Slv){
        _slv=Slv;
        _iterSlv=dynamic_cast<LinearSolverIter<Oper,V>*>(_slv.get());
        return *this;
    };

    OdeCrankNicolson<Oper,V> & withSolverTolerance(double Eps){_solverEps=Eps;return *this;};


    /// single Crank-Nicolson step
    V &step(V &Vec, double Tstart, double  Tstep){

        if(not _slv){
            _slv.reset(new LinearBiCGstab<LinSpaceMap<V>,V>(_one_hOp.get(),10));
            _slv->withTolerance(_solverEps==0.?1.e-12:_solverEps);
            _iterSlv=dynamic_cast<LinearSolverIter<Oper,V>*>(_slv.get());
        }

        derOde->update(Tstart+0.5*Tstep,z.get());
        *z=Vec;
        derOde->apply(0.5*Tstep,Vec,1,*z);

        _one_hOp->update(0.5*Tstep);
        _slv->compute(_one_hOp.get()); // update solver as needed

        // initial guess: 1+Tstep*Op*Vec, compute as -(-2*(1+0.5*Tstep*Op)*Vec + Vec)
        if(_iterSlv){
            Vec.axpy(-2.,*z);
            Vec.scale(-1.);
            Vec=_iterSlv->withInitialVector(&Vec).solve(*z);
        }
        else
            Vec=_slv->solve(*z);

        OdeStep<Oper,V>::nCallsStep=_iterSlv?_iterSlv->nApplys()+1:1;
        OdeStep<Oper,V>::nCalls+=OdeStep<Oper,V>::nCallsStep;
        return Vec;
    }

    /// return model vector argument of ODEstep
    const V & modelVec() const {return derOde->lhsVector();}

    unsigned int consistencyOrder() const {return 2;}
    double safetyFactor() const {return 0.5;}

    virtual std::string info() const {
        std::string s=OdeStep<Oper,V>::info();
        s+=(_iterSlv?", <iter>="+tools::str(double(_iterSlv->totIter())/_slv->nSolves(),3):"");
        return s;
    }

};

#endif // ODECRANKNICOLSON_H
