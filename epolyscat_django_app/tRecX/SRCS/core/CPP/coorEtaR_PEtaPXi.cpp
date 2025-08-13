#include "coorEtaR_PEtaPXi.h"
#include "abort.h"

std::vector<double> CoorEtaR_PEtaPXi::_fromRef(const std::vector<double> &PXiPEta) const{
    if(PXiPEta[0]==0. and PXiPEta[1]==0.)return {1.,0.};
    return {0.5*(PXiPEta[0]+PXiPEta[1]),(PXiPEta[0]-PXiPEta[1])/(PXiPEta[0]+PXiPEta[1])};
}

std::vector<double> CoorEtaR_PEtaPXi::_toRef(const std::vector<double> &EtaR) const{
    return {EtaR[1]*(1.+EtaR[0]),EtaR[1]*(1-EtaR[0])};
}

std::vector<double> CoorEtaR_PEtaPXi::_jacRefdCoor(const std::vector<double> &EtaR) const{
    if(EtaR[1]==0.)DEVABORT("no Jacobian at origin PXi,PEta=0,0");
    // dPXi/dEta, dPEta/dEta, dPXi/dR, dPEta/dR
    return {EtaR[1],-EtaR[1],1.+EtaR[0],1.-EtaR[0]};
}
