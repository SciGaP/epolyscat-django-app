// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "floquetAnalysis.h"

#include <sstream>
#include <iostream>
#include <string>

// for GetDiscSubregion
#include "operatorData.h"
//#include "operator.h" // this can go (minor changes)
#include "discretization.h"
#include "readInput.h"
#include "printOutput.h"
#include "eigenSolver.h"
#include "asciiFile.h"

using namespace std;

FloquetAnalysis::FloquetAnalysis(std::string FloquetHamiltonian,ReadInput & Inp)
    :floquetHamiltonian(FloquetHamiltonian)
{
    Inp.exclude("FloquetAnalysis","Eigen",ReadInput::anyName,"select");
    Inp.read("FloquetAnalysis","minE",eMin,"-2","(mutable) lower boundary for field free energies");
    Inp.read("FloquetAnalysis","maxE",eMax,"0","(mutable) upper boundary for field free energies");
    Inp.read("FloquetAnalysis","minOvr",minOvr,"0.1","(mutable) minimal overlap for a match");
    Inp.read("FloquetAnalysis","maxImag",maxImag,"1.e-1","(mutable) maximal modulus of imaginary part");
    Inp.read("FloquetAnalysis","maxErr",maxErr,"1.e-2","(mutable) maximal deviation between energies");
    int nFree;
    Inp.read("FloquetAnalysis","nFree",nFree,"10","(mutable) compare lowest nFree states of free problem");
    select="SmallReal["+tools::str(nFree)+"]";
    std::string dirRuns;
    Inp.read("FloquetAnalysis","compare",dirRuns,"","list of runs as in \"{Dir}[0001,0007,0009]\", Dir is omitted, use current top",1,"floquet");
    if(ReadInput::main.file()=="");
    if(dirRuns!=""){
        dir=dirRuns.substr(0,dirRuns.find("["))+"/";
        if(dir=="/"){
            dir=ReadInput::main.output();
            dir=dir.substr(0,dir.length()-ReadInput::inputCopy.length()-1);
        }
        runs=tools::splitString(tools::stringInBetween(dirRuns,"[","]"),',');
    }
}

void selectRange(double Emin, double Emax, double Imax,
                 const std::vector<std::vector<double>> Ener,
                 std::vector<std::pair<std::complex<double>, int>> & Sel){
    Sel.clear();
    for(size_t k=0;k<Ener[0].size();k++){
        if(Emin< Ener[0][k] and Ener[0][k] < Emax and std::abs(Ener[1][k])<Imax){
            Sel.push_back(std::pair<std::complex<double>, int>
                          (std::complex<double>(Ener[0][k],Ener[1][k]),k));
        }
    }
}

bool FloquetAnalysis::run(const Index* Idx)
{
    if(runs.size()==0)return false;

    ofstream floquet;
    floquet.open((ReadInput::main.output()+"floquet").c_str());
    floquet<<(Sstr+"# stable solutions from "+ReadInput::main.output()+","+dir+":"+runs)<<std::endl;
    floquet<<(Sstr+"# with eMin,eMax,minOvr,maxImag,maxErr: "+eMin+eMax+minOvr+maxImag+maxErr)<<std::endl;
    floquet<<"#    ReFloq  ImFloq  <Floq|Free> <Floq|Floq>   ReFree  "<<endl;

    // get the definition of a field free operator
    vector<string> terms = OperatorData::terms(floquetHamiltonian);
    string hamFieldFree = "";
    for(std::string t: terms){
        size_t idPos=t.find("<Id>");
        if(idPos<4)hamFieldFree+=t.substr(0,idPos)+t.substr(idPos+4);
    }
    if(hamFieldFree=="")ABORT("could not extract field-free from "+floquetHamiltonian);
    PrintOutput::message("The field free hamiltonian is "+hamFieldFree);

    // construct field free operator
    OperatorTree h0("HamFieldFree",hamFieldFree,Idx->child(0),Idx->child(0));
    EigenSolver free(eMin,eMax);
    free.withSelect(select);
    //    free.setMethod("Lapack");
    free.compute(&h0);
    free.eigenvalues();
    free.rightVectors();
    free.select(select);

    // select eigenvalues where all others have close values
    AsciiFile refEig(ReadInput::main.output()+"eig");
    std::vector<std::vector<double>> c;
    refEig.readCols(c);
    std::vector<std::pair<std::complex<double>,int>> eRef;
    selectRange(eMin,eMax,maxImag,c,eRef);

    // get eigenvectors matching selected eigenvalues
    for(std::string r: runs){
        AsciiFile compEig(dir+r+"/eig");
        c.clear();
        compEig.readCols(c);
        std::vector<std::pair<std::complex<double>,int>> eComp;
        selectRange(eMin,eMax,maxImag,c,eComp);
        for(int k=eRef.size()-1;k>=0;k--){
            size_t l=0;
            while(l<eComp.size() and abs(eRef[k].first-eComp[l].first)>maxErr)l++;
            if(l==eComp.size())eRef.erase(eRef.begin()+k);
        }
    }

    // loop through selected
    ifstream vecstream;
    vecstream.open((ReadInput::main.output()+"evec").c_str());
    Coefficients eVec(Idx);
    int nVec=0;
    for(auto e: eRef){
        while(eVec.read(vecstream,nVec==0) and e.second!=nVec++);
        if(not vecstream.good())DEVABORT("something wrong");

        // find block with largest norm
        int maxBlock(INT_MAX);
        double maxNrm=0;
        for(size_t k=0;k<Idx->childSize();k++){
            std::complex<double> nrm=Idx->child(k)->overlap()->matrixElement(*eVec.child(k),*eVec.child(k));
            if(maxNrm<std::abs(nrm)){
                maxNrm=std::abs(nrm);
                maxBlock=k;
            }
        }

        // find field-free state with largest overlap
        int nMatch(INT_MAX);
        double oMatch(0);
        for(size_t k=0;k<free.eigenvalues().size();k++){
            std::complex<double> ovr=Idx->child(k)->overlap()->matrixElement(*free.rightVectors()[k],*eVec.child(maxBlock));
            if(std::abs(ovr)>oMatch){
                oMatch=std::abs(ovr);
                nMatch=k;
            }
        }
        if(oMatch>minOvr)
            floquet<<e.first.real()<<", "<<e.first.imag()
                  <<", "<<std::norm(oMatch)/maxNrm<<", "<<maxNrm
                 <<", "<<free.eigenvalues()[nMatch].real()<<", "<<nMatch<<std::endl;
    }

    floquet.close();
    PrintOutput::message("selected Floquet energies on "+ReadInput::main.output()+"floquet",0,true);
    return true;
}
