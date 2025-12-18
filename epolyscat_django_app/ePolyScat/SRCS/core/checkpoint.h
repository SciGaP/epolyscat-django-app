// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include "tools.h"
#include "coefficients.h"
#include "threads.h"
#include "printOutput.h"

class Checkpoint{
    //    static constexpr std::string chptC="/chptC";
    std::string _dir;
    double _time;
public:
    Checkpoint(std::string Dir):_dir(Dir),_time(-DBL_MAX){
        std::ifstream ifil((Dir+"/chptC").c_str(),(std::ios_base::openmode) std::ios::beg|std::ios::binary);
        if(ifil.is_open())tools::read(ifil,_time);
    }
    bool operator()(){return _time>-DBL_MAX;}
    double time() const{return _time;}
    bool coefs(Coefficients * C){
        std::ifstream ifil((_dir+"/chptC").c_str(),(std::ios_base::openmode) std::ios::beg|std::ios::binary);
        if(ifil.is_open()){
            Coefficients* jC=Threads::join(*C);
            if(Threads::isMaster()){
                tools::read(ifil,_time);
                jC->read(ifil,false);
            }
            Threads::scatter(jC,*C);
            return true;
        }
        else
            return false;
    }
    void remove(){
        PrintOutput::message("remove checkpoint"+operator()());
        if(operator()())std::remove((_dir+"/chptC").c_str());
    }

    static void write(std::string WriteDir, double Time, const Coefficients* C){
        std::ofstream chptFile((WriteDir+"/chptC").c_str(),(std::ios_base::openmode) std::ios::beg);
        tools::write(chptFile,Time);
        C->write(chptFile,false);
        chptFile.close();
    }
};

#endif // CHECKPOINT_H
