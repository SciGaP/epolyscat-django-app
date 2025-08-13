#ifndef COORETARN_H
#define COORETARN_H

#include "coorSystem.h"

/// \ingroup Coordinates
/// \brief Eta.R coordinates with PXiPEta as reference
class CoorEtaR_PEtaPXi : public CoorSystem
{
protected:
    std::vector<double> _toRef(const std::vector<double> & EtaR) const;
    std::vector<double> _fromRef(const std::vector<double> & PXiPEta) const;
    std::vector<double> _jacRefdCoor(const std::vector<double> & EtaR) const;
public:
    CoorEtaR_PEtaPXi(std::string Name=""):CoorSystem("Eta.R","PXi.PEta",Name){}
};

#endif // COORETARN_H
