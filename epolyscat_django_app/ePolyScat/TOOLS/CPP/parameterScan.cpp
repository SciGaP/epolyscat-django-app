// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "parameterScan.h"

#include "mpiWrapper.h"
#include "readInput.h"
#include "parameters.h"
#include "printOutput.h"
using namespace std;

ParameterScan::ParameterScan(ReadInput &Inp)
{
    /// read multi-parameter ranges from input
    read(Inp,_names);
    outDir=Inp.outputTopDir();
}

void ParameterScan::scan(){

    bool first=true;
    vector<double> par,val;
    AsciiFile file(outDir+"scan");
    PrintOutput::title("scanning parameters: "+tools::str(_names));
    while(next(par)){
        valForPar(_names,par,val);

        if(first){
            first=false;
            vector<string> comm(1,"");
            for(int k=0;k<_names.size();k++)comm.back()+="  "+_names[k];
            for(int k=0;k<val.size();k++)comm.back()+="   val["+tools::str(k)+"]";
            file.writeComments(comm);
        }
        if(MPIwrapper::isMaster()){
            val.insert(val.begin(),par.begin(),par.end());
            file.writeRow(val,12);
        }
    }
    // print message, let it be repeated with the terminal prints for convenience (Final=true)
    PrintOutput::message("parameter scan on file "+file.name(),0,true);
}
