#ifndef BASISMATMO_H
#define BASISMATMO_H

#include <string>
#include "basisMatMatrix.h"

class BasisMO;
class BasisMatMO : public BasisMatMatrix
{
public:
    BasisMatMO(std::string Op, const BasisMO* IBas, const BasisMO*JBAS);
};

#endif // BASISMATMO_H
