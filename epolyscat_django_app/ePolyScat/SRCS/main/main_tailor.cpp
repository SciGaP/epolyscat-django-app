// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "readInput.h"
#include "printOutput.h"
#include "units.h"
#include "mpiWrapper.h"
#include "fieldTailor.h"

#include <dlib/optimization.h>

using namespace std;

int main(int argc, char* argv[]) {
    MPIwrapper::Init(argc,argv);

    //============================================================================================
    // input and part of setup
    //==============================================================================================

    ReadInput::openMain("",argc,argv);
    Units::setDefault("au");        // general default units
    ReadInput::main.setUnits("au"); // convert input to these units

    if(not MPIwrapper::isMaster())PrintOutput::off();

    // read the pulse range
    FieldTailor f(ReadInput::main);
    f.print();

    ReadInput::main.finish();

    // run all trajectories
    f.run();
    f.writeTrajectories();
}
