// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "time.h"
#include "platformSpecific.h"

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include <unistd.h>

namespace platformSpecific {

void current_utc_time(struct timespec *ts){

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts->tv_sec = mts.tv_sec;
    ts->tv_nsec = mts.tv_nsec;
#else
    //HACK for now
    clock_gettime(CLOCK_REALTIME, ts);
#endif

}

std::string current_host(){
#ifndef _STATIC_
    char cHost[255];
    if(gethostname(cHost,255)!=0)return "(unknown)";
    return std::string(cHost);
#else
    return "(unknown)";
#endif
}

}
