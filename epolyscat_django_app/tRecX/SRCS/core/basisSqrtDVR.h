#ifndef BASISSQRTDVR_H
#define BASISSQRTDVR_H

#include "basisIntegrableSqrt.h"
#include "basisDvr.h"

/// DVR basis: |Q|^(p/2) basDVR[n](Q) (limted to p=0,1 at present)
///
/// when an element boundary = 0, a Radau rule is used such that the boundary point does not appear in the rule
///<br> for use use with parabolic coordinates at m-quantum number >0 with condition p=m%2
class BasisSqrtDVR : public BasisDVR
{
    int _powSqrt;
public:
    virtual ~BasisSqrtDVR(){}
    BasisSqrtDVR(const std::string Def);

    /// values and derivatives at points X, ZeroOutside: set = 0, where X outside range of basis function
    void valDer(const std::vector<std::complex<double> > & X,
                std::vector<std::complex<double> > & Val,
                std::vector<std::complex<double> > & Der, bool ZeroOutside=false) const;
};

#endif // BASISSQRTDVR_H
