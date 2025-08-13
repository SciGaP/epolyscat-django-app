// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "newtonInterpolator.h"

#include "qtEigenDense.h"
#include "tools.h"
#include "readInput.h"
#include "wavefunction.h"
#include "coefficients.h"

using namespace std;
using namespace tools;
#include "eigenNames.h"


NewtonInterpolator::NewtonInterpolator(vector<complex<double> > &coordinates, vector<complex<double> > &fctValues) : supports(coordinates) {
    if(coordinates.size() != fctValues.size()) ABORT("Stop in NewtonInterpolator::NewtonInterpolator number of coordinates and function values does not agree!");
    // dividierte differenzen schema
    coeffs.push_back(fctValues[0]);
    vector<complex<double> > divDiffs = fctValues;
    while(divDiffs.size()>1) {
        for(unsigned int l=0; l<divDiffs.size()-1; l++) {
            divDiffs[l] = (divDiffs[l+1]-divDiffs[l])/(coordinates[l+1+coordinates.size()-divDiffs.size()] - coordinates[l]);
        }
        divDiffs.pop_back();
        coeffs.push_back(divDiffs[0]);
    }
}

vector<complex<double> >& NewtonInterpolator::getVals(vector<complex<double> > &coordinates, vector<complex<double> > &fctValues) const {
    // horner schema for evaluating the polynomial
    fctValues.assign(coordinates.size(),coeffs.back());
    for(int l = (int) supports.size()-2; l>-1; l--) {
        Map<ArrayXcd>(fctValues.data(),fctValues.size()) *= Map<ArrayXcd>(coordinates.data(),coordinates.size()) - supports[l];
        Map<ArrayXcd>(fctValues.data(),fctValues.size()) += coeffs[l];
    }
    return fctValues;
}

#include "wavefunction.h"
#include "discretization.h"

NewtonInterpolatorWF::NewtonInterpolatorWF(const Discretization *disc, int numberOfSupportPoints){

    // the expansion coefficients of the unique polynomial in newton-basis
    polynomialCoefficients.resize(numberOfSupportPoints);
    for(unsigned int i=0;i<polynomialCoefficients.size();i++)
        polynomialCoefficients[i] = new Wavefunction(disc);

    // the dividierte differenzen
    divDiffs.resize(numberOfSupportPoints);
    for(unsigned int i=0;i<divDiffs.size();i++)
        divDiffs[i] = new Wavefunction(disc);

    // the supportpoints
    times.resize(numberOfSupportPoints);

}

NewtonInterpolatorWF::NewtonInterpolatorWF(const Index* Idx, int numberOfSupportPoints){

    // the expansion coefficients of the unique polynomial in newton-basis
    polynomialCoefficients.resize(numberOfSupportPoints);
    for(unsigned int i=0;i<polynomialCoefficients.size();i++)
        polynomialCoefficients[i] = new Wavefunction(Idx);

    // the dividierte differenzen
    divDiffs.resize(numberOfSupportPoints);
    for(unsigned int i=0;i<divDiffs.size();i++)
        divDiffs[i] = new Wavefunction(Idx);

    // the supportpoints
    times.resize(numberOfSupportPoints);

}

void NewtonInterpolatorWF::computePolynomialCoefficients(const std::vector<Wavefunction *> &knownWfs){
    // dividierte differenzen schema
    // knownWfs->time represent the support points
    // knownWfs->coefs represent the known function values at the support points
    // the interpolation at a given time returns a wavefunction

    // some basic consistency checks
    if(knownWfs.size() != times.size()) ABORT("number of wave functions does not match number of times");
    for(unsigned int i=0;i<knownWfs.size()-1;i++) // check if times are monotonically increasing
        if(knownWfs[i]->time >= knownWfs[i+1]->time){
            cout << i << " " << knownWfs[i]->time << " " << knownWfs[i+1]->time << endl;
            ABORT("times are not in increasing order");
        }

    for(unsigned int i=0;i<divDiffs.size();i++){
        *(divDiffs[i]) = *(knownWfs[i]);
        times[i] = knownWfs[i]->time;
    }

    *(polynomialCoefficients[0]) = *(knownWfs[0]);

    int divDiffsCounter = divDiffs.size();
    int polCoeffCounter = 1;
    while(divDiffsCounter>1) {
        for(int l=0; l<divDiffsCounter-1; l++) {
            *(divDiffs[l]) -= *(divDiffs[l+1]);
            *(divDiffs[l]) *= ( 1.0 / ( knownWfs[l]->time - knownWfs[l+1+knownWfs.size()-divDiffsCounter]->time ) );
            //divDiffs[l] = (divDiffs[l+1]-divDiffs[l])/(coordinates[l+1+coordinates.size()-divDiffs.size()] - coordinates[l]);
        }
        divDiffsCounter--;
        *(polynomialCoefficients[polCoeffCounter]) = *(divDiffs[0]);
        polCoeffCounter++;
    }
}

NewtonInterpolatorWF::~NewtonInterpolatorWF(){
    for(unsigned int i=0;i<polynomialCoefficients.size();i++)
        delete polynomialCoefficients[i];
    for(unsigned int i=0;i<divDiffs.size();i++)
        delete divDiffs[i];
}

bool NewtonInterpolatorWF::inInterval(double Time){
    if(Time >= times[0]  and  Time <= times.back())return true;
    double eps=(abs(times[0])+abs(times.back()))*1.e-12;
    return times[0]<Time+eps and Time-eps<times.back();
}

void NewtonInterpolatorWF::getInterpolatedWF(double time, Wavefunction *result){
    // horner scheme for evaluating the polynomial
    if(not inInterval(time))ABORT("must not extrapolate: "+tools::str(time)
              +" not in ["+tools::str(times[0])+","+tools::str(times.back())+"]");

    *result = *(polynomialCoefficients.back());
    result->time=time;
    for(int l = (int) polynomialCoefficients.size()-2; l>-1; l--) {
        *result *= ( time - times[l] );
        *result += *(polynomialCoefficients[l]);
        result->time=time;
    }
}

void NewtonInterpolatorWF::timeDerivative(double time, Wavefunction *result){
    // use Horner scheme for evaluating derivative
    if(not inInterval(time))ABORT("must not extrapolate: "+tools::str(time)
              +" not in ["+tools::str(times[0])+","+tools::str(times.back())+"]");

    *result = *(polynomialCoefficients.back());
    result->coefs->scale(double(polynomialCoefficients.size()-1));
    result->time=time;
    // derivative: omit a_0 iterationd and
    // multiply all polynial coefficients by corresponding power bevor derivation
    for(int l = (int) polynomialCoefficients.size()-2; l>0; l--) {
        *result *= ( time - times[l] );
//        *result += *(polynomialCoefficients[l]);
        result->coefs->axpy(double(l),polynomialCoefficients[l]->coefs);
        result->time=time;
    }
}
