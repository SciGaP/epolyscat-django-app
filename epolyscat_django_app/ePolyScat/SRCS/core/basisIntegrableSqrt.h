#ifndef BASISINTEGRABLESQRT_H
#define BASISINTEGRABLESQRT_H

#include "basisIntegrable.h"

class BasisIntegrableSqrt : public BasisIntegrable
{
    int _powSqrt;
    const BasisIntegrable* _basInt;
public:
    static std::string inputToStrDefinition(BasisSetDef Def);
    static int powerSqrt(std::string Def);

    /// determines power in square root depending on position in Index tree, e.g interprete sqrt[m{Phi}%2]
    static void resolvePower(BasisSetDef &Def, std::vector<unsigned int> Branch, const std::vector<const Index*> Path);

    virtual ~BasisIntegrableSqrt(){}
    BasisIntegrableSqrt(const std::string Def);

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
