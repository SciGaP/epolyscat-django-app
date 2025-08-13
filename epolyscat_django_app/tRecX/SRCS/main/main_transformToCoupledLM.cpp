// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "tools.h"
#include <glob.h>
#include <string>
#include <vector>
#include "timer.h"

#ifdef _USE_BOOST_
#include <boost/algorithm/string/predicate.hpp>
#endif

#include "readInput.h"
#include "printOutput.h"

#include "discretization.h"
#include "coefficients.h"

#include "mpiWrapper.h"

#include "discretizationCoupledLM.h"
#include "coefficientsWriter.h"

using namespace std;

// Source: https://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
std::vector<std::string> glob(const std::string& pat){
    using namespace std;
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}

TIMER(all,)

int main(int argc, char* argv[]) {


    MPIwrapper::Init(argc,argv);

    START(all);

   // test for spectral parameter file


    ReadInput::openMain(string(argv[1])+"/inpc",argc,argv,false);
    ReadInput::main.setUnits("au");
    Units::setDefault("au");

    if(argc<2){
        PrintOutput::paragraph();
        PrintOutput::message("usage: >TransformToCoupledLM runDir/0123 [flags]");
        MPIwrapper::Finalize();
        exit(0);
    }
    if(argv[1][0]=='-')ABORT("first command line argument must be data directory, e.g. myRun/0137, found \""+string(argv[1])+"\"");

    PrintOutput::set(ReadInput::main.output()+"anaout");

    Discretization * D = Discretization::factory(ReadInput::main);
    D->print();

    PrintOutput::message("Setting up discretization (l1, l2, L, M)");
    DiscretizationCoupledLM disc(D);

    CoefficientsWriter base("base", "Rn1.Rn2");
    CoefficientsWriter transformed("transformed", "Rn1.Rn2");

    for(std::string path: glob(ReadInput::main.output()+"coeff/*")){
        std::string name = path.substr(path.find("coeff/")+6);
        double time;
        try{
            time = std::stod(name);
        }catch(std::exception& ex){
            // e.g. coeff/desc
            continue;
        }

        Coefficients c(D->idx());
        Coefficients cTransformed(disc.idx());

        std::ifstream inp(path, std::ios::binary);
        c.read(inp, false);
        inp.close();

        disc.mapFromParent()->apply(1., c, 0., cTransformed);
        
        base.write(time, c);
        transformed.write(time, cTransformed);

        PrintOutput::message("Transformed coeff@"+std::to_string(time));
    }

    STOP(all);

    PrintOutput::timerWrite();
    PrintOutput::title("done - "+string(argv[1]));
} // main
