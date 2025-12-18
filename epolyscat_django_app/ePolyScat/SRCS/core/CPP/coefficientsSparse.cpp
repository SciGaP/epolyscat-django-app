#include "../coefficientsSparse.h"
#include "index.h"

CoefficientsNull::CoefficientsNull():Coefficients(){}

CoefficientsSparse::CoefficientsSparse(const Coefficients & C):Coefficients()
{
    _cIndex=C.idx();
    if(C.isLeaf() and C.data()){
        if(not C.isZero()){
            reset(_cIndex);
            for(size_t k=0;k<_cIndex->size();k++)data()[k]=C.data()[k];
        }
    }
    for(size_t k=0;k<C.childSize();k++){
        if(not C.isZero())childAdd(new CoefficientsSparse(*C.child(k)));
        else              childAdd(new CoefficientsNull());
    }
}

Coefficients CoefficientsSparse::full() const {
    Coefficients res(idx());
    res=*this;
    return res;
}
