// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <functional>
#include "tools.h"
#include "stdio.h"
#include "qtAlglib.h"

namespace Integrate {

void test();

class Tools
{
protected:
    std::string kind,kindInf;
public:
    std::vector<std::vector<unsigned int> >nQuad;
    const double epsRel;
    const double epsAbs;
    const unsigned int maxLevel;
    unsigned int peakLevel;
    unsigned int funcCount;

public:
    std::string message; ///< prepend to error message

    /// integrator parameters
    std::string str() const {return tools::str(epsRel,2)+"/"+tools::str(epsAbs)+" ["+tools::str(nQuad[0])+"] "+kind+"/"+kindInf;}


    static const unsigned int defaultNQuad;

    static std::map<std::string,std::map<int,std::vector<std::vector<double> > > > rules;

    /// \brief Integrate parameters for integration
    /// \param AccRel    desired relative accuracy
    /// \param AccAbs    desired absolute accuracy
    /// \param NQuad     number of quadrature points
    /// \param Kind      quadrature rule on finite intervals
    /// \param KindInf   quadrature rule on infinite intervals
    Tools( double AccRel=1.e-12, double AccAbs=1.e-12,
           std::vector<std::vector<unsigned int> > NQuad=std::vector<std::vector<unsigned int> >(0),
           std::string Kind="GaussLegendre",
           std::string KindInf="GaussLaguerre",
           unsigned int MaxLevel=20)
        : epsRel(AccRel),epsAbs(AccAbs),nQuad(NQuad),kind(Kind),kindInf(KindInf),maxLevel(MaxLevel),
          peakLevel(0),funcCount(0),message("--- no message defined for integration ---"){}

    /// \brief quadRule quadrature points and weights for kind,kindInf of integrator
    /// \param A        lower boundary
    /// \param B        upper boundary
    /// \param QPoin    points, size=order
    /// \param QWeig    weights
    void quadRule(const std::vector<double> &Vol, std::vector<double> & QPoin,std::vector<double> & QWeig) const;

    /// add a new quadrature rule for [0,1] or [0,infty)
    void addQuadRule(const std::string Name, const std::vector<double> & QPoin, const std::vector<double> & QWeig);

};

// =====================================================================================================================
// function definition templates need to be included here, so functions can be created and compiled as needed
// function definition full specialization MUST NOT be included here, will cause "multiple definitions"

using namespace std;


void showVol(std::vector<std::vector<double> > Vol,std::string Text="Volume");///< boundaries of integration volume
double measureVol(vector<vector<double> > Vol); ///< size of integration volume

vector<vector<double> > & NextSubvolume(const std::vector<std::vector<double> > &Vol, std::vector<std::vector<double> > &Sub);

template <class ReturnType,class ArgType,class IntType>
ReturnType NDim(IntType& Int,
                const vector<vector<double> > & Vol,
                const std::function<ReturnType(const std::vector<ArgType> & )> Func,
                const std::vector<ArgType> & Params=std::vector<ArgType>(0),
                const std::vector<ArgType> & XVals=std::vector<ArgType>(0))
{
    if(Vol.size()==XVals.size()){
        Int.funcCount++;
        return Func(XVals);
    }

    unsigned int nQ=10;//Integrate::defaultNQuad;
    if(Int.nQuad.size()!=0)nQ=Int.nQuad[XVals.size()][0];
    vector<double> qPoin(nQ);
    vector<double> qWeig(qPoin.size());
    Int.quadRule(Vol[XVals.size()],qPoin,qWeig);

    vector<ArgType>xvals(XVals);
    xvals.push_back(qPoin[0]);
    ReturnType result(NDim(Int,Vol,Func,Params,xvals)*qWeig[0]);
    for(unsigned int k=1;k<qPoin.size();k++){
        xvals.back()=qPoin[k];
        result+=NDim(Int,Vol,Func,Params,xvals)*qWeig[k];
    }
    return result;
}

template <class ReturnType> double maxAbsVal(const ReturnType & X);
template <class ReturnType,class ArgType,class IntType>
ReturnType Recursive( IntType &I,
                      const vector<std::vector<double> > Vol,
                      const function<ReturnType(const vector<ArgType> &)> Func,
                      const vector<ArgType> Params=std::vector<ArgType>(0),
                      unsigned int Level=0, ReturnType Previous=ReturnType())
{

    if(Level==0){
        I.funcCount=0;
        Previous = I.nDim(Vol,Func,Params);
    }

    ReturnType integral(Previous);

    // systematically run through sub-volumes
    integral*=0.;
    vector<vector<double> > sub;
    vector<ReturnType> subInts;
    while (NextSubvolume(Vol,sub).size()>0){
        subInts.push_back(I.nDim(sub,Func,Params));
        integral+=subInts.back();
    }

    // check for relative and absolute errors
    if(maxAbsVal<ReturnType>(integral-Previous)<max(I.epsAbs,I.epsRel*maxAbsVal<ReturnType>(integral)))return integral;

    if(Level>I.maxLevel){
        cout<<"\nFAILURE in integrate\n "+I.message<<endl;
        showVol(Vol);
        cout<<maxAbsVal<ReturnType>(integral)<<" - "<<maxAbsVal<ReturnType>(Previous)<<" >? "<<I.epsAbs<<endl;
        if(measureVol(Vol)*1.e5<maxAbsVal(integral))cout<<"!!! integrand may be singular !!!"<<endl;
        ABORT("exceeded maximal recursion depth="+tools::str(I.maxLevel)
              +"\ncheck whether integrand has singularities or non-analyticities"
                +"\nrelax accuracy or increase order: "+I.str());
    }
    // get accurate integrals on each subvolume
    integral*=0.;
    for(unsigned int k=0;k<subInts.size();k++){
        // adjust tolerance according to importance of sub-volume
        integral+=Recursive(I,NextSubvolume(Vol,sub),Func,Params,Level+1,subInts[k]);

    }
    return integral;
}

template <class ReturnType,class ArgType>
class Int:public Integrate::Tools
{
public:
    Int( double AccRel=1.e-12, double AccAbs=1.e-12,
         std::vector<std::vector<unsigned int> > NQuad=std::vector<std::vector<unsigned int> >(0),
         std::string Kind="GaussLegendre",
         std::string KindInf="GaussLaguerre")
        : Integrate::Tools(AccRel,AccAbs,NQuad,Kind,KindInf){}

    ReturnType nDim(const std::vector<std::vector<double > > Vol,
                    const std::function<ReturnType(const std::vector<ArgType> &)> Func,
                    const std::vector<ArgType> Params=std::vector<ArgType>(0)
            ) {return Integrate::NDim<ReturnType,ArgType,Int>(*this,Vol,Func,Params);}

    ReturnType recursive(const std::vector<std::vector<double> > Vol,
                         const std::function<ReturnType(const std::vector<ArgType> &)> Func,
                         const std::vector<ArgType> Params=std::vector<ArgType>(0)
            ) {return Integrate::Recursive<ReturnType,ArgType,Int>(*this,Vol,Func,Params);}
};

}


#endif // INTEGRATE_H
