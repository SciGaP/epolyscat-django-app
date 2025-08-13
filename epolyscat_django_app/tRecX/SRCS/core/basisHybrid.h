#ifndef BASISHYBRID_H
#define BASISHYBRID_H

#include "basisVector.h"

/// wrapper for BasisVector use on hybrid axes
class BasisHybrid : public BasisVector
{
public:
    BasisHybrid():BasisVector(0){_name="Hybrid";}
    BasisHybrid(int Size):BasisVector(Size){}
    std::string strDefinition() const override {return "Hybrid:"+tools::str(size());}

};

#endif // BASISHYBRID_H
