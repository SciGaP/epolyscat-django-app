#include "coorParabolic.h"
#include "abort.h"

std::vector<double> CoorParabolic::_fromRef(const std::vector<double> &EtaR) const{
    return {EtaR[1]*(1.+EtaR[0]),EtaR[1]*(1-EtaR[0])};
}

std::vector<double> CoorParabolic::_toRef(const std::vector<double> &PXiPEta) const{
    if(PXiPEta[0]==0. and PXiPEta[1]==0.)return {1.,0.};
    return {(PXiPEta[0]-PXiPEta[1])/(PXiPEta[0]+PXiPEta[1]),(PXiPEta[0]+PXiPEta[1])/2};
}

std::vector<double> CoorParabolic::_jacRefdCoor(const std::vector<double> &PXiPEta) const{
    if(PXiPEta[0]==0. and PXiPEta[1]==0.)DEVABORT("no Jacobian at origin PXi,PEta=0,0");
    double den=2/std::pow(PXiPEta[0]+PXiPEta[1],2);
    return {PXiPEta[1]*den,0.5,-PXiPEta[0]*den,0.5};
}
