// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "tools.h"

#include <iostream>
#include <fstream>
#include "mpiWrapper.h"
#include "spectrumLinCom.h"
#include "spectrumPlot.h"
#include "readInput.h"
#include "units.h"
#include "math.h"
#include "printOutput.h"

using std::ofstream;
using std::ios_base;

/// Linearly combine amplitudes for amplitudes at different phiCEO
///
/// applicable in the perturbative limit
///
/// Usage, spectrum at specific phiCEO:
/// <br>
/// LinCom VinayData/0_0 -plot=total -phiCEO=2*pi/3 -ampl=VinayData/0_pi2 -channel=0 -interpolate
///
/// Usage, determine phi with maximal yield from yield in [kmin,kmax]=[0.72,0.79]:
/// <br>
/// LinCom VinayData/0_0 -plot=total -phiCEO=maxPhi[0.72,0.79] -ampl=VinayData/0_pi2 -channel=0 -interpolate -lambda=800
///
/// Usage with tRecX generated amplitudes
/// <br>
/// LinCom 71scr/0001 -ampl=2 -phiCEO=maxPhi
///
/// for 71scr/0001 and .../0002
/// k-range will be determined from spectra as peak between two highest peaks
/// all pulse parameters will be taken from original 71scr/0001 ../0002 inputs
int main(int argc, char* argv[]) {
    MPIwrapper::Init(argc,argv);
    Units::setDefault("au");        // general default units
    ReadInput::main.setUnits("au"); // convert input to these units
    ReadInput::openMain("",argc,argv);
    ReadInput::main.suppressLinp(); // do NOT write a linp file for this executable (keep run's linp)
    PrintOutput::set(ReadInput::main.outputTopDir()+"/outLinCom");

    double lambda;
    ReadInput::main.read("Laser","lambda(nm)",lambda,"0","wave length",1,"lambda");
    if(lambda==0)ABORT("give wave length in nm, e.g. -lambda=800");
    Units::addUnit("OptCyc|time",lambda*1.e-9/physics::speed_of_light);

    SpectrumPlot plt(ReadInput::main.outputTopDir(),"total",true,{});
    SpectrumLinCom lc;
    lc.plotInputSpectra(plt);
    lc.plotAverage(plt);
    PrintOutput::paragraph();
    PrintOutput::message(Str("delays for lambda=","")+lambda+" nm, (OptCyc="+Units::convert(1.,"OptCyc")+" au)");
    PrintOutput::title("done");
//    ReadInput::main.finish(); // original input is used, finish won't work as most items are not addressed
    ReadInput::main.writeDoc();
}
