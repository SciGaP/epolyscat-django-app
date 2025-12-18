#include "basisSqrtDVR.h"

#include "basisIntegrableSqrt.h"

BasisSqrtDVR::BasisSqrtDVR(const std::string Def)
    :BasisDVR(Def.substr(Def.find("*")+1))
{
    _powSqrt=BasisIntegrableSqrt::powerSqrt(Def);
    _name=Def.substr(0,Def.find("*")+1)+_name;
}

void BasisSqrtDVR::valDer(const std::vector<std::complex<double> > &X,
                          std::vector<std::complex<double> > &Val,
                          std::vector<std::complex<double> > &Der, bool ZeroOutside) const
{
    BasisDVR::valDer(X,Val,Der,ZeroOutside);
    if(_powSqrt==0)return;
    if(_powSqrt>1)DEVABORT("power of square root > 1 not implmented: "+str(0));

    std::vector<double> sqrt,qSqrt2;
    for(auto x: X){
       sqrt.push_back(std::sqrt(std::abs(x)));
       qSqrt2.push_back(sqrt.back()==0.?0:0.5/sqrt.back());
    }

    for(size_t n=0,kn=0;n<size();n++){
        // we keep normalizations same as underlying BasisDVR
        double scalN=_dvrX[_nBeg+n]==0.?0.:1./std::sqrt(std::abs(_dvrX[_nBeg+n]));
        for(size_t k=0;k<X.size();k++,kn++){
            Der[kn]=(sqrt[k]*Der[kn]+qSqrt2[k]*Val[kn])*scalN;
            Val[kn]=(sqrt[k]*Val[kn])*scalN;
        }
    }
}
