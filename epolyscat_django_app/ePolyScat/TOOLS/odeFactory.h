#ifndef ODEFACTORY_H
#define ODEFACTORY_H

#include "odeStep.h"
#include "odeRK4.h"
#include "odeIntegral.h"
#include "odeCrankNicolson.h"
#include "odeImplicitRK.h"
#include "tools.h"

template<class Der,class V>
OdeStep<Der,V>* odeFactory(std::string Method, Der* D, double Eps=0.){
    OdeStep<Der,V>* res=0;
    if(Method=="RK4")res=new OdeRK4<Der,V>(D);
    else if(Method=="euler")res=new OdeIntegral<Der,V>(D);
    else if(Method=="CrankNicolson")res=new OdeCrankNicolson<Der,V>(D,Eps*0.01);
    else if(Method.find("RK[")==0)res=new OdeImplicitRK<Der,V>(D,std::max(1e-11,std::min(1.e-7,Eps*0.01)),
                                                               tools::stringInBetween(Method,"RK[","]",true));

    if(not res)ABORT("unknown ODE step method "+Method+", available: RK4,euler,CrankNicolson");
    return res;
}

#endif // ODEFACTORY_H
