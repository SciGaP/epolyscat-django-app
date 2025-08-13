// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "abort.h"
#include "../TOOLS/toolsHeader.h"
#include "mpiWrapper.h"
#include "readInput.h"
#include <map>
#include <iostream>
void abortFile(const std::string file,int line, const std::string mess, int signal) {
    std::cerr <<file.substr(file.rfind("/")+1)<<":"<<line;
    if(mess!=""){
        if(MPIwrapper::Size()>1)
            std::cerr<<"\n<"<<MPIwrapper::Rank()<<"> "<<mess<<"\n";
        else
            std::cerr<<"\n"<<mess << "\n";
    }
    std::cerr<< std::endl << std::flush;
#ifdef _DEVELOP_
    signal=1;
#endif
    if(MPIwrapper::Size()>1)
        MPIwrapper::Abort(signal);
    else
        abort();
}
void devAbortFile(const std::string file,int line, const std::string mess) {
#ifdef _DEVELOP_
    abortFile(file,line,"(Developer) "+mess,1);
#else
    abortFile(file,line,"(Developer) "+mess+"\n --- Developer assistance may be required --- ",0);
#endif
}


std::map<std::string,unsigned int> abortCnt;
void abortCountDown(const std::string file,int line,const std::string mess, unsigned int Cnt) {
    if(abortCnt.count(file+mess)==0)abortCnt[file+mess]=Cnt;
    if (abortCnt[file+mess]==0)abortFile(file,line,mess,1);
    abortCnt[file+mess]--;
}


/// \brief activate by flag -DEBUGnew[=true]
bool newCode(){return ReadInput::main.flag("DEBUGnew","use new code during changes");}
/// \brief activate by flag -DEBUGold[=true]
bool oldCode(){return ReadInput::main.flag("DEBUGold","use old code during changes");}
bool oldEva(){return ReadInput::main.flag("DEBUGeva","use old evaluator during changes");}


/** @} */
