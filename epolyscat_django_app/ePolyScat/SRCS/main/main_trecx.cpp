// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "tRecX.h"
#include "farm.h"
#include "inputGenerator.h"
#include "threads.h"
#include "printOutput.h"

#include "readInput.h"
#include "problemTrecX.h"
#include "log.h"

int main(int argc, char* argv[]) {
    //    Sstr+"try this A"+boost::math::sph_bessel(1,0.27557017407092*20)+Sendl;

    std::srand(time(NULL)); // seed for (likely) random upon repeated runs

    MPIwrapper::Init(argc,argv);

    // try get from problem data base
    std::shared_ptr<ProblemTrecX> prob(new ProblemTrecX(argc,argv));

    Farm farm;
    if(prob->nInps()==1)farm=Farm(prob);
    else                farm=Farm(InputGenerator::factory(argc,argv));

    if(MPIwrapper::isMaster())PrintOutput::set("cout");
    //    ProblemTree prob(farm.inp());

    // split into (possibly multi-process) threads
    farm.fork();

    // loop through runs in present thread
    while(farm.next()){
        if( not farm.runCompleted()){
            tRecX::run_trecx(farm.inp(),argc,argv);
            tools::sleepSeconds(0.1);
        }
    }
    farm.join();
    MPIwrapper::Barrier();
    tools::sleepSeconds(0.1);
    farm.status();

    MPIwrapper::setCommunicator(MPIwrapper::worldCommunicator());
    if(MPIwrapper::isMaster(MPIwrapper::worldCommunicator())){
        PrintOutput::warningList();
        std::cout<<"\n === tRecX done =================================================="<<std::endl;
    }
    MPIwrapper::Barrier();
    MPIwrapper::Finalize();
    exit(0);
}
