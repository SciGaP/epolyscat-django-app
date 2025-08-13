#ifndef COEFFICIENTSSPARSE_H
#define COEFFICIENTSSPARSE_H

#include "coefficients.h"

class CoefficientsSparse : public Coefficients
{
public:
    virtual ~CoefficientsSparse(){};
    CoefficientsSparse(const Coefficients & C);
    Coefficients full() const;
};

class CoefficientsNull: public Coefficients
{
public:
    virtual ~CoefficientsNull(){};
    bool notNull() const {return false;}
    CoefficientsNull();//:Coefficients(){}
};
#endif // COEFFICIENTSSPARSE_H
