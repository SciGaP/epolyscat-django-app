#ifndef ALGEBRAPIECWISE_H
#define ALGEBRAPIECWISE_H

#include "algebra.h"
#include <memory>

/// sum of algebras on disjoint intervals (-infty,q0,...qN,infty)
/// where q0,...,qN are the InternalBoundaries
class AlgebraPiecewise : public Algebra
{
protected:
    std::vector<double> _qInternal;
    std::vector<std::shared_ptr<Algebra>> _algs;
    std::complex<double> _scale;
    AlgebraPiecewise():_scale(1.){}
public:
    /// compose by Algebra definitions Def, separated by InternalBoundaries
    AlgebraPiecewise(std::vector<double> & InternalBoundaries, std::vector<std::string> Defs);

    /// set scaling factor applied to all pieces
    void withScale(std::complex<double> Scale){_scale=Scale;}

    std::complex<double> val(std::complex<double> Q) const;
    std::complex<double> integral(std::complex<double> Q0, std::complex<double> Q1) const;
    std::vector<std::complex<double>> nonAnalyticQ() const;
    std::string str() const;
};

#endif // ALGEBRAPIECWISE_H
