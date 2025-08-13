#ifndef BASISDYSON_H
#define BASISDYSON_H

#include "basisAbstract.h"
#include "vectorValuedFunction.h"
#include "basisMO.h"

class BasisDyson :  public BasisAbstract, public VectorValuedFunction
{
    std::string _source;
    const VectorValuedFunction* _mo;
    std::vector<std::vector<double>> _dysonC;

    BasisDyson(const VectorValuedFunction* MO, const std::vector<std::vector<double>> & DysC, std::string Source);
public:
    static BasisDyson* read(const VectorValuedFunction *MO, std::string DysonCoefs);
    std::vector<std::complex<double>> operator()(std::vector<double> X) const;
    std::string  coordinates() const {return _mo->coordinates();}
    unsigned int length() const {return _mo->length();}
    unsigned int size() const {return length();}
    std::string source() const {return _source;}
    void print() const;
    bool isOrthonormal() const {return true;}
    std::string strDefinition() const{return "undefined";}

};

#endif // BASISDYSON_H
