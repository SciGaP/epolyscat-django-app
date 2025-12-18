// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "scanEigenvalue.h"

#include "readInput.h"
#include "printOutput.h"
#include "derivativeFlat.h"
//#include "operator.h"
#include "parameters.h"

using namespace std;

vector<complex<double> > ScanEigenvalue::resolv;
vector<vector<complex<double> > >ScanEigenvalue::eval;
vector<vector<Coefficients*> >ScanEigenvalue::evec;
EigenSubspace * ScanEigenvalue::eigSub=0;

ScanEigenvalue::~ScanEigenvalue(){
    for(int k=0;k<evec.size();k++)
        for(int l=0;l<evec[k].size();l++)delete evec[k][l];
}


ScanEigenvalue::ScanEigenvalue(ReadInput &Inp):ParameterScan(Inp)
{
    if(not Inp.found("ScanEigenvalues"))return;

    // add parameters to resetable table (defaults to value=1.)
    for(int k=0;k<_names.size();k++)Parameters::addResetable(_names[k]);

    while(not Inp.endCategory("ScanEigenvalues",resolv.size()+1)){
        double reGuess,imGuess,reVal,imVal;
        Inp.read("ScanEigenvalues","reGuess",reGuess,ReadInput::noDefault,"real part(s) of guess target energy",resolv.size()+1);
        Inp.read("ScanEigenvalues","imGuess",imGuess,"0","imaginary part of guess target energy",resolv.size()+1);
        Inp.read("ScanEigenvalues","reResolv",reVal,tools::str(reGuess),"real part of E0 in preconditioner (H0-E0)^-1",resolv.size()+1);
        Inp.read("ScanEigenvalues","imResolv",imVal,tools::str(imGuess),"real part of E0 in preconditioner (H0-E0)^-1",resolv.size()+1);
        resolv.push_back(complex<double>(reGuess,imGuess));
        eval.push_back(vector<complex<double> >(1,complex<double>(reVal,imVal)));
    }
    if(resolv.size()==0)ABORT("need at least one reGuess");

    evec.resize(eval.size());
    eigSub=new EigenSubspace(Inp);
    Parameters::setSpecial();
    setFunction(eigenvaluesAtPar);
}

void ScanEigenvalue::print() const {

    MultiParam::print();

    PrintOutput::title("Eigenvalues by subspace iteration");
    //    PrintOutput::lineItem("size",int(e0.size()));
    PrintOutput::newLine();
    eigSub->print();
    PrintOutput::newLine();

    PrintOutput::newRow();
    PrintOutput::rowItem("reE0");
    PrintOutput::rowItem("imE0");
    PrintOutput::rowItem("reGuess");
    PrintOutput::rowItem("imGuess...");
    for(int k=0;k<resolv.size();k++){
        PrintOutput::newRow();
        PrintOutput::rowItem(resolv[k].real());
        PrintOutput::rowItem(resolv[k].imag());
        for(int l=0;l<eval[k].size();l++){
            PrintOutput::rowItem(eval[k][l].real());
            PrintOutput::rowItem(eval[k][l].imag());
        }
    }
    PrintOutput::paragraph();

}

void ScanEigenvalue::eigenvaluesAtPar(const std::vector<std::string> &ParName, const std::vector<double> &ParVal, std::vector<double> &Result)
{
    for(int k=0;k<ParName.size();k++)Parameters::reset(ParName[k],ParVal[k]);
    Parameters::updateSpecial();
    Result.clear();
    for(int k=0;k<eval.size();k++){
        eigSub->eigen(vector<complex<double> >(1,resolv[k]),eval[k],evec[k]);
        Result.push_back(eval[k][0].real());
        Result.push_back(eval[k][0].imag());
    }
}
