// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "integrate.h"
#include "useMatrix.h"
#include "orthopol.h"
#include "vectorReal.h"
#include "qtEigenDense.h"

//#include "basisFunction.h"

using namespace std;

/// show the present integration volume
void Integrate::showVol(std::vector<std::vector<double> > Vol,std::string Text) {
    std::cout<<Text+": ";
    for(unsigned int k=0;k<Vol.size();k++)std::cout<<"["<<Vol[k][0]<<","<<Vol[k][1]<<"] ";
    std::cout<<std::endl;
}

double Integrate::measureVol(vector<vector<double> > Vol){
    double m=1;
    for(unsigned int k=0;k<Vol.size();k++)m*=(Vol[k][1]-Vol[k][0]);
    return m;
}

/// initialize to first if Sub=empty on entry
/// return Sub=empty vector if cannot increment
vector<vector<double> > & Integrate::NextSubvolume(const std::vector<std::vector<double> > &Vol, std::vector<std::vector<double> > &Sub){

    // for infinite intervals: shift by vol[k][2]*fac
    const double fac=0.5;

    // reached emtpy volume, return empty
    if(Vol.size()==0){
        Sub.clear();
        return Sub;
    }

    // enter with empty sub-volume, initialize
    if(Sub.size()==0){
        Sub=Vol;
        for(unsigned int k=0;k<Vol.size();k++){
            if(Vol[k][1]>DBL_MAX/2.)
                Sub[k][1]=Vol[k][0]+Vol[k][2]*fac;
            else if(Vol[k][0]<-DBL_MAX/2.)
                Sub[k][1]=Vol[k][1]+Vol[k][2]*fac;
            else
                Sub[k][1]=0.5*(Vol[k][0]+Vol[k][1]);
        }
        return Sub;
    }

    // try increment last
    if(Sub.back()[0]==Vol.back()[0]){
        Sub.back()[1]=Vol.back()[1];
        if(Vol.back()[1]>DBL_MAX/2.)
            Sub.back()[0]=Vol.back()[0]+Vol.back()[2]*fac;
        else if(Vol.back()[0]<-DBL_MAX/2.)
            Sub.back()[0]=Vol.back()[1]+Vol.back()[2]*fac;
        else
            Sub.back()[0]=0.5*(Vol.back()[0]+Vol.back()[1]);
        return Sub;
    }

    // create leading term
    vector<vector<double> > vol(Vol);
    Sub.pop_back();
    vol.pop_back();

    // could not increment, return empty
    if(NextSubvolume(vol,Sub).size()==0)return Sub;

    // successful increment, reset trailing interval
    Sub.push_back(Vol.back());
    if(Vol.back()[1]>DBL_MAX/2.)
        Sub.back()[1]=Vol.back()[0]+Vol.back()[2]*fac;
    else if(Vol.back()[0]<-DBL_MAX/2.)
        Sub.back()[1]=Vol.back()[1]+Vol.back()[2]*fac;
    else
        Sub.back()[1]=0.5*(Vol.back()[0]+Vol.back()[1]);
    return Sub;

}


namespace Integrate{
template<>double maxAbsVal<double>(const double & X ){return abs(X);}
template<>double maxAbsVal<complex<double> >(const complex<double> & X ){return abs(X);}
template<>double maxAbsVal<UseMatrix>(const UseMatrix & Mat){return Mat.maxAbsVal();}
template<>double maxAbsVal<VectorReal>(const VectorReal & Vec){return Vec.maxAbsVal();}
template<>double maxAbsVal<Eigen::MatrixXcd>(const Eigen::MatrixXcd & Mat){return Mat.lpNorm<Eigen::Infinity>();}
}

/// \brief IntExample - how to use IntTools and templates IntNDim, IntRecursive
///
class IntExample:public Integrate::Tools
{

public:
    /// \brief Integrate parameters for integration
    /// \param AccRel    desired relative accuracy
    /// \param AccAbs    desired absolute accuracy
    /// \param NQuad     number of quadrature points
    /// \param Kind      quadrature rule on finite intervals
    /// \param KindInf   quadrature rule on infinite intervals
    IntExample( double AccRel=1.e-12, double AccAbs=1.e-12,
                std::vector<std::vector<unsigned int> > NQuad=std::vector<std::vector<unsigned int> >(0),
                std::string Kind="GaussLegendre",
                std::string KindInf="GaussLaguerre")
        : Integrate::Tools(AccRel,AccAbs,NQuad,Kind,KindInf){}

    double nDim(const std::vector<std::vector<double> > Vol,
                const std::function<double(const std::vector<double> &)> Func,
                const std::vector<double> Params=std::vector<double>(0)
            ) {return Integrate::NDim<double,double,IntExample>(*this,Vol,Func,Params);}

    double recursive(const std::vector<std::vector<double> > Vol,
                     const std::function<double(const std::vector<double> &)> Func,
                     const std::vector<double> Params=std::vector<double>(0)
            ) {return Integrate::Recursive<double,double,IntExample>(*this,Vol,Func,Params);}
};

const unsigned int Integrate::Tools::defaultNQuad=8;
map<string,map<int,vector<vector<double> > > >Integrate::Tools::rules;

void Integrate::Tools::addQuadRule(const string Name, const std::vector<double> &QPoin, const std::vector<double> &QWeig){
    if(QPoin.size()!=QWeig.size())ABORT("number of points and weight differ");
    if(rules[Name][QPoin.size()].size()==0){
        rules[Name][QPoin.size()]=vector<vector<double> >(1);
        rules[Name][QPoin.size()].push_back(QPoin);
        rules[Name][QPoin.size()].push_back(QWeig);
    }
}

void Integrate::Tools::quadRule(const std::vector<double>& Interval, std::vector<double> &QPoin, std::vector<double> &QWeig) const{
    // get points and weights for quadrature for interval
    if(QPoin.size()<1 )ABORT("need at least 1 quadrature point");
//    if(QPoin.size()>99)ABORT("use at most 99 quadrature points");
    if(Interval[1]> DBL_MAX/2. and Interval[2]<=0.)ABORT("specify positive scale for infinite upper boundary: "+tools::str(Interval));
    if(Interval[0]<-DBL_MAX/2. and Interval[2]>=0.)ABORT("specify negative scale for infinite lower boundary: "+tools::str(Interval));

    double scal(1.);
    string thisKind;
    if(abs(Interval[0])<DBL_MAX/2. and abs(Interval[1])<DBL_MAX/2.){
        thisKind=kind;
        scal=Interval[1]-Interval[0];
    }
    else if(Interval[0]<-DBL_MAX/2. or Interval[1]>DBL_MAX/2.){
        thisKind=kindInf;
        if(Interval.size()<3)ABORT("infinite interval needs scaling as third element in interval, is: "+tools::str(Interval,8,","));
        scal=Interval[2];
    }
    else
        ABORT("cannot integrate over interval (-infty,infty): "+tools::str(Interval));

    vector<vector<double> > Rule(rules[thisKind][QPoin.size()]);
    vector<double> xVec,wVec,wFun;
    double shift=0.,scale=1.;
    if(Rule.size()==0){
        // need to put rule into table
        Rule.resize(2);
        if     (thisKind=="lobatto")
        {
            OrthogonalLegendre().quadratureWithEnds(QPoin.size(),xVec,wVec);
            shift=1.;
            scale=0.5;
        }
        else if(thisKind=="radau")
        {
            OrthogonalLaguerre().quadratureWithEnds(QPoin.size(),xVec,wVec);
            for(unsigned int k=0;k<xVec.size();k++)wFun.push_back(OrthogonalLaguerre().weight(xVec[k]));
        }
        else if(thisKind=="GaussLegendre")
        {
            OrthogonalLegendre().quadrature(QPoin.size(),xVec,wVec);
            shift=1.;
            scale=0.5;
        }
        else if(thisKind=="GaussLaguerre"){
            OrthogonalLaguerre().quadrature(QPoin.size(),xVec,wVec);
            for(unsigned int k=0;k<xVec.size();k++)wFun.push_back(OrthogonalLaguerre().weight(xVec[k]));
        }

        else if(thisKind=="index"){
            // basis quadrature always for [0,1]
            for(unsigned int k=0;k<QPoin.size();k++){
                scale=1./scal;
                shift=0.;
                xVec.push_back(double(k));
                wVec.push_back(1.);
            }
        }
        else
            ABORT("quadrature not implemented: "+thisKind);

        // divide by weight function (if any)
        for(unsigned int k=0;k<wFun.size();k++)wVec[k]/=wFun[k];

        UseMatrix x(xVec.size(),1),w(wVec.size(),1);
        for(unsigned int k=0;k<xVec.size();k++){
            x(k)=(xVec[k]+shift)*scale;
            w(k)=wVec[k]*scale;
        }

        for(size_t i = 0; i<x.size(); i++) {
            Rule[0].push_back(x(i).real());
            Rule[1].push_back(w(i).real());
        }
        rules[thisKind][QPoin.size()]=Rule;
    }

    // check whether correct length rule was returned
    if(QPoin.size()!=Rule[0].size())
        ABORT("incorrect Quadratur rulen length: "+tools::str(QPoin.size())+" vs. "+tools::str(Rule[0].size()));

    // scale to present interval
    for (unsigned int k=0;k<QPoin.size();k++){
        if    (Interval[0]>-DBL_MAX/2)QPoin[k]=Interval[0]+scal*Rule[0][k];
        else if(Interval[1]<DBL_MAX/2)QPoin[k]=Interval[1]+scal*Rule[0][k];
        else ABORT("cannot handle interval: "+tools::str(Interval));
        QWeig[k]=  abs(scal)*Rule[1][k];
    }
    return;
}

double Func1(const vector<double>&X)
{
    double expon=0;
    for(unsigned int k=0;k<X.size();k++)expon-=X[k]*double(k+1);
    return exp(expon);
}
void Integrate::test(){
    IntExample I1;

    vector<vector<double> >vol(1,vector<double>(2,0.));
    vol[0][1]=1.;
    cout<<"Int[0,1] exp(-x)"<<setprecision(14)<<endl;
    cout<<"    fixed "<<I1.nDim(vol,Func1)<<endl;
    cout<<"recursive "<<I1.recursive(vol,Func1)<<endl;
    cout<<"    exact "<<1.-exp(-1.)<<endl;
    cout<<I1.str()<<endl;

    vol.push_back(vol[0]);
    cout<<"Int[0,1][0,1] exp(-x-2y)"<<setprecision(14)<<endl;
    cout<<"    fixed "<<I1.nDim(vol,Func1)<<endl;
    cout<<"recursive "<<I1.recursive(vol,Func1)<<endl;
    cout<<"    exact "<<(1.-exp(-1.))*0.5*(1.-exp(-2.))<<endl;
    cout<<" "<<I1.str()<<endl;

    vol.push_back(vol[0]);
    vol.front()[1]=3.;
    cout<<"Int[0,3][0,1][0,1] exp(-x-2y-3z)"<<setprecision(14)<<endl;
    cout<<"    fixed "<<I1.nDim(vol,Func1)<<endl;
    cout<<"recursive "<<I1.recursive(vol,Func1)<<endl;
    cout<<"    exact "<<(1.-exp(-3.))*(1.-exp(-2.))/2.*(1.-exp(-3.))/3.<<endl;
    cout<<" "<<I1.str()<<endl;

    vol.front()[1]=DBL_MAX;
    vol.front().push_back(1.);
    vol.pop_back();
    vol.pop_back();
    IntExample I2(1.e-12,1.e-12);
    I2.nQuad.push_back(vector<unsigned int>(2,4));
    cout<<"Int[0,infty][0,1] exp(-x)"<<setprecision(14)<<endl;
    cout<<"    fixed "<<I2.nDim(vol,Func1)<<endl;
    cout<<"recursive "<<I2.recursive(vol,Func1)<<endl;
    cout<<"    exact "<<1.<<endl;
    cout<<" "<<I2.str()<<endl;

}
