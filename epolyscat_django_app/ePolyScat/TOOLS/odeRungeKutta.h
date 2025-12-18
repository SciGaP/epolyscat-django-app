// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ODERUNGEKUTTA_H
#define ODERUNGEKUTTA_H

#include "odeStep.h"

#include <complex>
#include <vector>
#include <complex>
#include "abort.h"
#include "unistd.h"
#include "str.h"
#include "printOutput.h"

#include "butcherTableau.h"

class TestExp{
    std::complex<double> phas;
public:
    TestExp(double Phas=0.):phas(std::exp(std::complex<double>(0.,Phas))){}
    void apply(const std::complex<double> A, const std::complex<double> X,
               const std::complex<double> B, std::complex<double> & Y){
        Y*=B;
        Y+=A*X*phas;
    }
    void update(double Time){} //
    std::complex<double> exact(std::complex<double> Y0, double Step){
        return Y0*exp(phas*Step);
    }
    std::complex<double> rhsVector(){return 0.;}


};

template<class Der,class V>
class OdeRungeKutta: public OdeStep<Der,V>{
    using OdeStep<Der,V>::derOde; // C++? surprising it does not know where it is derived from
    std::vector<std::vector<double> > a;
    std::vector<double> b,c;
    std::vector<V> vecK;
public:
    OdeRungeKutta(Der*D,const ButcherTableau & BuTab)
        :OdeStep<Der,V>(BuTab.name(),D),
          a(BuTab.a()),b(BuTab.b()),c(BuTab.c())
    {
        OdeStep<Der,V>::_consistency=BuTab.consistency();
        // check whether method is explicit
        for(int i=0;i<a.size();i++)
            for (int j=i;j<a.size();j++)
                if(a[i][j]!=0.)ABORT(Str("only for explicit Runge Kutta methods, row=")+i+j+a[i]);
        OdeStep<Der,V>::nCallsStep=a.size();
        vecK.resize(b.size(),D->rhsVector());
    }

    /// explicit RungeKutta step
    V &step(V &Vec, double Tstart, double  Tstep){

        vecK.back()=Vec;
        for(int i=0;i<a.size();i++){
            vecK[i]=vecK.back();
            for(int j=0;j<i;j++){
                if(a[i][j]!=0.)vecK[i]+=Tstep*a[i][j]*vecK[j];
            }

            derOde->update(Tstart+Tstep*c[i]);      // set time to t0+h*c[i]
            derOde->apply(1.,vecK[i],0.,vecK[i]);
            Vec+=Tstep*b[i]*vecK[i];
        }

        OdeStep<Der,V>::nCalls+=nApplyStep();
        return Vec;
    }

    /// return model vector argument of ODEstep
    const V & modelVec() const {return vecK[0];}

    unsigned int consistencyOrder() const {return OdeStep<Der,V>::_consistency;}
    unsigned int nApplyStep() const {return a.size();}

    double numConsistency(){
        double h=1.;
        double e0=1.,e1=1.;
        while (e0+e1>1.e-8 and h>1.e-14){
            e0=e1;
            std::complex<double> y(1.);
            step(y,0,h);
            e1=abs(y-derOde->exact(1.,h));
            h*=0.5;
        }
        if(h<=1.e-14)ABORT("could not reach errror < 1.e-8");
        return log(e0/e1)/log(2.)-1;
    }

    double stability(){
        // binary search for stability
        double h0=0.,h1=1.;
        std::complex<double> y1=1.;
        step(y1,0.,h1);
        while(abs(y1)<1.){
            h1*=2.;
            y1=1.;
            step(y1,0.,h1);
        }
        while(h1-h0>1.e-14){
            y1=1.;
            step(y1,0.,(h0+h1)*0.5);
            if(abs(y1)<1.)h0=(h0+h1)*0.5;
            else          h1=(h0+h1)*0.5;
        }
        return h0;
    }

    static bool test(){
        std::vector<std::string> meth={"RK4","Butcher67"};
        for(int k=0;k<meth.size();k++){
            for(int l=0;l<3;l++){
                TestExp der(3.1415926*(-2-l)*0.25);
                OdeRungeKutta<TestExp,std::complex<double> >ode(&der,ButcherTableau(meth[k]));
                Str("Arg(degree), method, order(theo,num):")+45*(-2-l)+ode.name()
                        +ode.consistencyOrder()+ode.numConsistency()
                        +"stability:"+ode.stability()
                        +Str::print;
            }
        }
        return true;
    }
};

#endif // ODERUNGEKUTTA_H
