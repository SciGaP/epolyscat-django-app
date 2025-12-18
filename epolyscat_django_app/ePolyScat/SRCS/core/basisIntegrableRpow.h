#ifndef BASISINTEGRABLERPOW_H
#define BASISINTEGRABLERPOW_H

#include "basisIntegrable.h"

/// @brief behaves as R^l at R=0
///
/// on interval [0,R1]:  Q^l*polynomial (i.e. Jacobi polynomial with B=2*l)
///<br> on other interval:       polynomial
///<br> Basis order and size reduced by l, minimal order 3
///<br> value of last function at b(Q=R1)=1
class BasisIntegrableRpow: public BasisIntegrable
{
    int _powQ;
    const BasisIntegrable* _basInt; // ownership is handled by BasisAbstract::factory
    Eigen::MatrixXd _trans;

public:
    static std::string inputToStrDefinition(BasisSetDef Def);

    /// determines power of basis depending on position in Index tree, e.g interprete Rpow[l{Eta}]
    ///
    /// R^l is only multiplied in interval starting from 0
    ///<br> else the basis after the Rpow*-factor is employed
    ///<br> if RpowL* is specified, power will be derived from unique Eta-axis
    static void resolvePower(BasisSetDef &Def,std::vector<unsigned int> Branch, const std::vector<const Index*> Path);

    virtual ~BasisIntegrableRpow(){}
    BasisIntegrableRpow(const std::string Def);
    std::string strDefinition() const;

    unsigned int size()  const {return _basInt->size();}
    unsigned int order() const {return _basInt->order();};
    void quadRule(int N, std::vector<double> &QuadX, std::vector<double> &QuadW) const
    {_basInt->quadRule(N,QuadX,QuadW);};

    /// values and derivatives at points X, ZeroOutside: set = 0, where X outside range of basis function
    void valDer(const std::vector<std::complex<double> > & X,
                std::vector<std::complex<double> > & Val,
                std::vector<std::complex<double> > & Der, bool ZeroOutside=false) const;

};

#endif // BASISINTEGRABLESQRT_H
