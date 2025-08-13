// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include <chrono>
#include <ctime> // Probably better with recent C++ functionality, but StackOverflow loves C headers
#include <iomanip>

#include "log.h"
#include "mpiWrapper.h"
#include "readInput.h"
#include "memInfo.h"

Log Log::main;

void Log::off(){
    enabled = false;
}

void Log::on(){
    enabled = true;
}

void Log::push(std::string name, std::string file, int line, std::string func){
    if(!enabled) return;

    stack.push_back(name);
    log("DEBUG", file, line, func) << "PUSH" << std::flush;
    mem_info(file, line, func, true) << "PUSH";
}

void Log::pop(std::string file, int line, std::string func){
    if(!enabled) return;

    if(stack.size() == 0) return;
    stack.pop_back();
    log("DEBUG", file, line, func) << "POP" << std::flush;
    mem_info(file, line, func, true) << "POP";
}

void Log::flush(){
    if(!enabled) return;

    if(output.is_open()){
        output << std::flush;
        output.close();
    }

    output.open(ReadInput::main.output()+"log_"+std::to_string(MPIwrapper::Rank())+".csv",
            initial ? std::ofstream::out : std::ofstream::app);
}

std::ostream& Log::log(std::string severity, std::string file, int line, std::string func){
    if(!enabled) return output;

    if(not output.is_open()) flush();
    else output << std::endl;

    if(initial){
        output <<"time, timestamp, severity, stack, location, function, message" << std::endl;
        initial = false;
    }

    long timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
            ).count();

//    auto t = std::time(0);

    // HACK - my compiler does not know put_time...?
#ifdef _LOG_
//    auto tm = *std::localtime(&t);
//    output << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
#endif
    output << "," <<  timestamp << "," << severity << ",";

    for(size_t i=0; i<stack.size(); i++){
        output << stack[i];
        if(i < stack.size() - 1) output << "-";
    }

    output << "," << file.substr(file.rfind("/")+1) << ":" << line << "," << func << ",";

    output << std::flush;
    return output;
}

std::ofstream& Log::log(){
    if(!enabled) return output;

    output << std::flush;
    return output;
}
    
std::ofstream& Log::mem_info(std::string file, int line, std::string func, bool force){
    if(!enabled) return output;

    std::map<std::string, double> mem_info = getMemInfo();

    if(!force &&
            std::abs(last_rss-mem_info["current_rss"]) < LOG_MEM_STEP &&
            std::abs(last_peak_rss-mem_info["peak_rss"]) < LOG_MEM_STEP)
        return dummy;

    last_rss = mem_info["current_rss"];
    last_peak_rss = mem_info["peak_rss"];

    log("MEM_INFO", file, line, func) << std::setprecision(11) << "free(GiB): " <<
        mem_info["free_ram"] << "/" << mem_info["total_ram"] << " --- RSS(GiB): " <<
        mem_info["current_rss"] << "/" << mem_info["peak_rss"] << " --- ";
    return output;
}
