#ifndef LINEARBICGSTAB_H
#define LINEARBICGSTAB_H

#include "LinearSolverIter.h"
#include <memory>
#include <complex>
#include "abort.h"
#include "timeCritical.h"
#include "str.h"


/// BiCGstab solver for linear system Op x = b, with restart
template<class Oper, class V>
class LinearBiCGstab: public LinearSolverIter<Oper,V>{
    using LinearSolver<Oper,V>::_op;
    using LinearSolverIter<Oper,V>::_maxIter;
    using LinearSolverIter<Oper,V>::_nIter;
    using LinearSolverIter<Oper,V>::_nApplys;
    using LinearSolverIter<Oper,V>::_eps;
    using LinearSolverIter<Oper,V>::_guess;
    using LinearSolverIter<Oper,V>::_residues;

    static constexpr double defaultEpsRestart=1.e-6;
    double _epsRestart;

    std::shared_ptr<V> r0,rj,xj,sj,pj,Asj,Apj;
public:
    LinearBiCGstab(const Oper *Op,int MaxIter /** maximal number of iterations */)
        :LinearSolverIter<Oper,V>(Op,MaxIter),_epsRestart(defaultEpsRestart)
    {
        _residues.assign(_maxIter,0.);
    };

    /// set threshold for resatsrt
    LinearBiCGstab& withRestartThreshold(double EpsRestart){_epsRestart=EpsRestart; return *this;}


    void compute(const Oper *Op){
        // already set up
        if(r0)return;
        timeCritical::suspend();
        r0.reset(new V(Op->lhsVector()));
        rj.reset(new V(*r0));
        xj.reset(new V(*r0));
        sj.reset(new V(*r0));
        pj.reset(new V(*r0));
        Asj.reset(new V(*r0));
        Apj.reset(new V(*r0));
        timeCritical::resume();
    }

    V& solve(const V& B) const
    {
        if(not r0)DEVABORT("call LinearBiCGstab::compute(...) before solve(...)");

        LinearSolverIter<Oper,V>::_nSolves++;
        //-----------------------------------------------------
        // algorithm as in tsurff.pdf (bicgstab.tex)
        //------------------------------------------------------
        // 1: rj=b; A.apply(-1,xj,1,rj)
        _nApplys=0;
        *rj=B;
        if(_guess){
            *xj=*_guess;
            _op->apply(-1,*xj,1,*rj);
            _nApplys++;
        }
        else
            *xj*=0.;

        std::complex<double> rjr0(0.); // [4] rjr0=0
        // 3:for j<maxIter:
        for(_nIter=0;_nIter<_maxIter;_nIter++){
            LinearSolverIter<Oper,V>::_totIter++;

            // 18:    if rjr0 small:
            if(std::abs(rjr0)<_epsRestart){
                *r0=*rj;// 19:r0=rj
                // 20: pj=rj; rjr0=rj.r0
                *pj=*rj;
                rjr0=rj->dotProduct(*r0);
            }

            // 4:    A.apply(1,pj,0,Apj); aj=rjr0/Apj.r0
            _op->apply(1,*pj,0,*Apj);
            _nApplys++;
            std::complex<double> aj=rjr0/Apj->dotProduct(*r0);

            // 5:    sj=rj;  sj.axpy(-aj,Apj)
            *sj=*rj;
            sj->axpy(-aj,*Apj);

            // 6: if(sj small)
            double sjNrmsq=std::abs(sj->dotProduct(*sj));
            if(sjNrmsq<_eps*_eps){
                _residues[_nIter]+=sjNrmsq;
                xj->axpy(aj,*pj); // 7: xj.axpy(aj,pj)
                break; // 8:
            }

            _op->apply(1,*sj,0,*Asj); //[10] A.apply(1,sj,0,Asj)
            _nApplys++;
            // 10: oj=Asj.sj/Asj.Asj
            std::complex<double> oj=Asj->dotProduct(*sj)/Asj->dotProduct(*Asj);
            // 11: xj.axpy(aj,pj); xj.axpy(oj,sj)
            xj->axpy(aj,*pj);
            xj->axpy(oj,*sj);
            // 12: rj=sj; rj.axpy(-oj,Asj)
            *rj=*sj;
            rj->axpy(-oj,*Asj);

            // 13:    if(rj small) break
            double rjNrmsq=std::abs(rj->dotProduct(*rj));
            _residues[_nIter]+=rjNrmsq;
            if(rjNrmsq<_eps*_eps)break;

            // 16:    bj=aj/(oj*rjr0); rjr0=rj.r0; bj*=rjr0;
            std::complex<double> bj=aj/(oj*rjr0);
            rjr0=rj->dotProduct(*r0);
            bj*=rjr0;
            // 17:    sj=pj; A.apply(-bj*oj,sj,bj,pj); pj+=rj
            *sj=*pj;
            _op->apply(-bj*oj,*sj,bj,*pj);
            _nApplys++;
            *pj+=*rj;
        }
        return *xj;
    }
};



#endif // LINEARBICGSTAB_H
