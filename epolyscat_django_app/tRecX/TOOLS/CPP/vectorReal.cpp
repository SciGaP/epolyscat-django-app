// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "vectorReal.h"
#include "tools.h"

using namespace std;

VectorReal VectorReal::operator-(const VectorReal & Other) const {
    VectorReal c(*this);
    double* a=c.data();
    const double * b=Other.data();
    for(;a!=c.data()+size();a++,b++)*a-=*b;
    return c;
}

VectorReal & VectorReal::operator+=(const VectorReal & Other) {
    double* a=data();
    const double * b=Other.data();
    for(;a!=data()+size();a++,b++)*a+=*b;
    return *this;
}
VectorReal & VectorReal::operator-=(const VectorReal & Other) {
    double* a=data();
    const double * b=Other.data();
    for(;a!=data()+size();a++,b++)*a-=*b;
    return *this;
}

VectorReal & VectorReal::operator*=(double A){
    for(double* a=data();a!=data()+size();a++)*a*=A;
    return *this;
}
VectorReal VectorReal::operator*(double A) const {
    VectorReal b(*this);
    for(double* a=b.data();a!=b.data()+size();a++)*a*=A;
    return b;
}

VectorReal & VectorReal::axpy(double A, const VectorReal & X){
    double* a=data();
    const double * b=X.data();
    for(;a!=data()+size();a++,b++)*a+=*b*A;
    return *this;
}

double VectorReal::maxAbsVal() const {
    double maxAbs=0;
    for(const double* a=data();a!=data()+size();a++)maxAbs=max(maxAbs,std::abs(*a));
    return maxAbs;
}
double VectorReal::dot(const VectorReal & B) const {
//    if(size()!=B.size())ABORT("vector length do not match");
    double res=0.;
    for(int k=0;k<size();k++)res+=data()[k]*B[k];
    return res;
}

double VectorReal::normSqu() const {
    double nrm=0.;
    for(const double* a=data();a!=data()+size();a++)nrm+=pow(*a,2);
    return nrm;
}

VectorReal & VectorReal::purge(double Eps){
    for(int k=0;k<size();k++)
        if(abs(data()[k])<maxAbsVal()*Eps and maxAbsVal()>Eps)data()[k]=0.;
    return *this;
}

void VectorReal::Test(){
    VectorReal a,b;
    a.assign(10,3.);
    b.resize(a.size(),2.);
    b+=a;
    cout<<tools::str(b,3,", ")<<endl;
    b*=2.;
    cout<<tools::str(b,3,", ")<<endl;
}
