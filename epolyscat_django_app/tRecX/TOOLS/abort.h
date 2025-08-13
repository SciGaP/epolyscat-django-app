// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef _ABORT_H_
#define _ABORT_H_
#include <map>
#include <iostream>

#define DEVABORT(message) {devAbortFile(__FILE__,__LINE__, message);exit(1);}
#define ABORT(message) {abortFile(__FILE__,__LINE__,message);exit(1);}
//#define DEVABORT(message)  exit(0)
//#define ABORT(message) exit(0)
#define COUNTDOWN(message,cnt) abortCountDown(__FILE__,__LINE__,message,cnt)


/** @defgroup Abort Run time
 *  @ingroup Tools
  * \brief controlled abort, timing
  * @{
  */
/// \file

/// \ingroup Abort
/// \fn static void abortFile(const std::string file,int line, const std::string mess="")
/// \brief universal abort function, best use through macro ABORT(message-string)
void abortFile(const std::string file,int line, const std::string mess="",int signal=0);
void devAbortFile(const std::string file,int line, const std::string mess="");

/// \ingroup Abort
/// \fn
/// \brief countdown abort function, best use through macro COUNTDOWN(message-string, counts)
void abortCountDown(const std::string file,int line,const std::string mess, unsigned int Cnt);


/// \brief activate by flag -DEBUGnew[=true]
bool newCode();
/// \brief activate by flag -DEBUGold[=true]
bool oldCode();
bool oldEva();

#endif

/** @} */
