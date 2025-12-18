#include "algebraPiecewise.h"

AlgebraPiecewise::AlgebraPiecewise(std::vector<double> & InternalBoundaries, std::vector<std::string> Defs)
    :AlgebraPiecewise()
{
    _qInternal=InternalBoundaries;
    if(not std::is_sorted(_qInternal.begin(),_qInternal.end()))ABORT(Sstr+"need interval boundries sorted, got: "+_qInternal)
            if(_qInternal.size()+1!=Defs.size())ABORT("need one more algebra then internal interval boundaries");
    for(auto def: Defs)_algs.push_back(std::shared_ptr<Algebra>(new Algebra(def)));
    for(auto alg: _algs)if(not alg->isAlgebra())ABORT("cannot interprete algebra string "+alg->definition());
}

std::complex<double> AlgebraPiecewise::val(const std::complex<double> Q) const{
    size_t n=0;
    while(n<_qInternal.size() and Q.real()>_qInternal[n])n++;
    return _scale*_algs[n]->val(Q);
}

std::vector<std::complex<double>> AlgebraPiecewise::nonAnalyticQ() const{
    std::vector<std::complex<double>> res;
    for (auto q: _qInternal)res.push_back(q);
    return res;
}

std::complex<double> AlgebraPiecewise::integral(const std::complex<double> Q0, const std::complex<double> Q1) const{
    std::complex<double> res(0.);
    if(_qInternal.front()>Q0.real())res+=_algs.front()->integral(Q0,_qInternal.front());
    for(size_t n=1;n<_qInternal.size();n++){
        double qlo=std::max(_qInternal[n-1],Q0.real());
        double qup=std::min(_qInternal[n  ],Q1.real());
        if(qlo<qup)res+=_algs[n]->integral(qlo,qup);
    }
    if(_qInternal.back()<Q1.real())res+=_algs.back()->integral(_qInternal.back(),Q1);
    return res;
}


std::string AlgebraPiecewise::str() const {
    std::string res="{"+_algs.front()->definition()+"}";
    for(size_t k=0;k<_qInternal.size();k++){
        res+=tools::str(_qInternal[k],3)+"{"+_algs[k+1]->definition()+"}";
    }
    return res;
}
