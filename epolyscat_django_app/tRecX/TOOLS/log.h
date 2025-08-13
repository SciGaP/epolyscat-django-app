// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef LOG_H
#define LOG_H

/**
 * \file
 * Simple macro based logging for debugging purposes. Every process logs to individual file: log_<rank>.csv.
 * Files are closed and opened at every macro (however, before the new content is written).
 *
 * File format (csv):
 * \code{.csv}
 * 2018-03-01 9:00:01, <unix timestamp>, DEBUG/INFO/..., <stack>, main.cpp:555, function, <message>
 * \endcode
 * where <stack> can be set using LOG_PUSH and LOG_POP (see below).
 * 
 * Basic usage:
 * \code{.cpp}
 * LOG_[DIWE](x<<" == "<<y<<"? "<<(x == y));
 * \endcode
 * 
 * where DIWE stands for DEBUG, INFO, WARN, ERROR.
 *
 * Continue a started log message:
 * \code{.cpp}
 * LOG_[DIWE]("Part1... ");
 * LOG("Part2: "<<x);
 * \endcode
 *
 * Stack usage:
 * \code{.cpp}
 * LOG_PUSH("func1");
 * 
 * LOG_D("abc");
 *
 * LOG_PUSH("func2");
 *
 * LOG_D("xyz");
 *
 * LOG_POP();
 * LOG_POP();
 * \endcode
 *
 * producing:
 * \code{.csv}
 * ..., ..., DEBUG, func1, ..., ..., abc
 * ..., ..., DEBUG, func1-func2, ..., ..., xyz
 * \endcode
 *
 * Mem info:
 * \code{.cpp}
 * LOG_MEM("additional...");
 *
 * // Only log if current or peak rss have increased by at least LOG_MEM_STEP 
 * LOG_MEM_IF_INCREASED("additional...");
 * \endcode
 *                                            
 * will log free/total ram and current/peak rss as a MEM_INFO message. LOG_MEM_IF_INCREASED is automatically called when
 * LOG_POP is called.
 *
 * Stack-based information on performace and memory consumption can be accessed using
 * SCRIPTS/log.py
 *
 */

#include <vector>
#include <fstream>

#define LOG_MEM_STEP 0.01

#ifdef _LOG_

/// Enable/Disable logging
#define LOG_ON() Log::main.on();
#define LOG_OFF() Log::main.off();

/// Push on top of stack
#define LOG_PUSH(name) {Log::main.push(name, __FILE__, __LINE__, __func__); Log::main.flush();}

/// Pop from top of stack
#define LOG_POP()      {Log::main.pop(__FILE__, __LINE__, __func__); Log::main.flush();}

/// Append to last log
#define LOG(content)   {Log::main.log()<<content; Log::main.flush();}

/// DEBUG level logging
#define LOG_D(content) {Log::main.log("DEBUG", __FILE__, __LINE__, __func__)<<content; Log::main.flush();}

/// INFO level logging
#define LOG_I(content) {Log::main.log("INFO",  __FILE__, __LINE__, __func__)<<content; Log::main.flush();}

/// WARN level logging
#define LOG_W(content) {Log::main.log("WARN",  __FILE__, __LINE__, __func__)<<content; Log::main.flush();}

/// ERROR level logging
#define LOG_E(content) {Log::main.log("ERROR", __FILE__, __LINE__, __func__)<<content; Log::main.flush();}

/// DEBUG level logging of current free/total ram and current/peak rss
#define LOG_MEM(content) {Log::main.mem_info(__FILE__, __LINE__, __func__, true)<<content; Log::main.flush();}

/// Only log MEM_INFO if current or peak rss have changed by at least LOG_MEM_STEP
#define LOG_MEM_IF_INCREASED(content) {Log::main.mem_info(__FILE__, __LINE__, __func__, false)<<content; Log::main.flush();}

#else

/// Enable/Disable logging
#define LOG_ON() ;
#define LOG_OFF() ;

#define LOG_PUSH(name) ;
#define LOG_POP()      ;
#define LOG(content)   ;
#define LOG_D(content) ;
#define LOG_I(content) ;
#define LOG_W(content) ;
#define LOG_E(content) ;

#define LOG_MEM(content) ;
#define LOG_MEM_IF_INCREASED(content) ;

#endif

#include <iostream>

class Log{
    bool enabled;
    bool initial;

    std::ofstream output;
    std::ofstream dummy;
    std::vector<std::string> stack;

    double last_rss;
    double last_peak_rss;

public:
    static Log main;

    Log():enabled(true), initial(true),last_rss(-LOG_MEM_STEP), last_peak_rss(-LOG_MEM_STEP){}

    void off();
    void on();
    void push(std::string name, std::string file, int line, std::string func);
    void pop(std::string file, int line, std::string func);
    void flush();
    std::ostream& log(std::string severity, std::string file, int line, std::string func);
    std::ofstream& log();
    std::ofstream& mem_info(std::string file, int line, std::string func, bool force);

};

#endif
