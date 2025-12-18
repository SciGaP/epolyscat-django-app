// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "coefficientsWriter.h"

#include <iomanip>

#include "folder.h"
#include "tools.h"
#include "readInput.h"
#include "mpiWrapper.h"
#include "averageOverAxes.h"
#include "wavefunction.h"
#include "index.h"
#include "timeCritical.h"

std::string CoefficientsWriter::default_name;
std::string CoefficientsWriter::default_averageAxes;
double CoefficientsWriter::default_store;

CoefficientsWriter* CoefficientsWriter::instance(){
    static CoefficientsWriter* inst = 0;
    if(inst==0){
        inst = new CoefficientsWriter(default_name, default_averageAxes, default_store);
    }
    return inst;
}

void CoefficientsWriter::read(ReadInput& Inp){
    Inp.read("CoefficientsWriter", "name", default_name, "", "output name of coefficients, this enables writing");
    Inp.read("CoefficientsWriter", "averageAxes", default_averageAxes, "", "Average over certain axes, e. g. \"Rn1.Rn2\"");
    Inp.read("CoefficientsWriter", "store",default_store,"0","minimal interval for file output");
}

CoefficientsWriter::CoefficientsWriter(std::string Name, std::string AverageAxes, double Store): 
    wfOutput(0),average(0),name(Name), averageAxes(AverageAxes), store(Store), last_store(-DBL_MAX/10.){

    if(name != ""){
        folder = ReadInput::main.output()+name+"/";

        if(MPIwrapper::isMaster()){
            if(not folder::exists(folder)){
                folder::create(folder);
            }
            if(not folder::exists(folder)){
                ABORT("Could not create folder "+folder);
            }
        }

    }
}

static void writeIndex(const Index* Idx, std::ofstream& Out){
    Out << std::setprecision(15);

    std::vector<std::string> hierarchy = tools::splitString(Idx->hierarchy(), '.');
    for(int i=0; i<hierarchy.size(); i++){
        Out << hierarchy[i];
        if(i<hierarchy.size()-1) Out << ",";
    }
    Out << std::endl;

    for(const Index* idx = Idx->firstLeaf(); idx!=0; idx=idx->nextLeaf()){
        std::vector<double> physical;
        for(const Index* i=idx; i->parent()!=0; i=i->parent()){
            physical.push_back(i->physical());
        }

        for(int k=physical.size()-1; k>=0; k--){
            Out << physical[k];
            if(k!=0) Out << ",";
        }
        Out << std::endl;
    }
}

void CoefficientsWriter::write(double Time, const Coefficients& Coeff){
    if(name=="") return; // Disabled
    if(Time - last_store < store) return;
    last_store = Time;

    if(not MPIwrapper::isMaster()){
        MPIwrapper::Barrier();
        return;
    }

    if(not wfOutput){
        timeCritical::suspend();
        const Index* idx = Coeff.idx();

        if(averageAxes!=""){
            average = new AverageOverAxes(Coeff.idx(), tools::splitString(averageAxes, '.'));
            idx = average->Idx;
            wfOutput = new Wavefunction(average->Idx);
        }
        else
            wfOutput = new Wavefunction(Coeff.idx());

        std::ofstream desc(folder+"desc.csv");
        writeIndex(idx, desc);

        timeCritical::resume();
    }
    
    std::ofstream output(folder+std::to_string(Time), std::ios::binary);

    if(average)average->average(Coeff, *wfOutput->coefs);
    else      *wfOutput->coefs=Coeff;
    wfOutput->time=Time;
    wfOutput->write(output, false);

    MPIwrapper::Barrier();
}

void CoefficientsWriter::disable(){
    name = "";
}
