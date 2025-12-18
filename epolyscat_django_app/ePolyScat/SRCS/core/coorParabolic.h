#ifndef COORPARABOLIC_H
#define COORPARABOLIC_H


#include "coorSystem.h"

/// \ingroup Coordinates
/// \brief parabolic coordinates
class CoorParabolic: public CoorSystem{
protected:
    std::vector<double> _toRef(const std::vector<double> & PXiPEta) const;
    std::vector<double> _fromRef(const std::vector<double> & EtaR) const;
    std::vector<double> _jacRefdCoor(const std::vector<double> & PXiPEta) const;
public:
    CoorParabolic(std::string Name=""):CoorSystem("PXi.PEta","Eta.R",Name){}
};

#endif // COORPARABOLIC_H
