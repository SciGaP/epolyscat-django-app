#include "runDir.h"
#include <algorithm>
#include <ios>
#include <sstream>
#include <cmath>        // std::pow
#include <iomanip> // std::setfill

#include "abort.h"
#include "folder.h"
#include "stringTools.h"

void RunDir::updateRoot(std::string Root){
    if(_root!=Root){
        _root=Root;
        if(_root.back()!='/')_root.push_back('/');
        updateDigits(digits());
    }
}

void RunDir::updateDigits(int Digits){
    _digits=Digits==0?digits():std::max(Digits,4);
}

std::string RunDir::dir() const {
    if(_cur>=std::pow(10,_digits))return "";
    std::stringstream ss;
    ss<<std::setfill('0') << std::setw(_digits) << _cur;
    return _root+ss.str()+"/";
}

RunDir &RunDir::nextUnused(){
    std::string res;
    // find lowest non-existing
    while(dir()!="" and folder::exists(dir()))_cur++;
    return *this;
}

int RunDir::digits(){
    // incomplete check: six or fewer-digits directores 0-10
    // once we switch to C++17, we will use std::filesystem instead
    int savRun=_cur;
    int savDigs=_digits;
    int res=4;
    for(_digits=7;_digits>4 and res==4;_digits--){
        for(_cur=0;_cur<10 and res==4;_cur++){
            if(folder::exists(dir())){
                res=_digits;
            }
        }
    }
    _digits=savDigs;
    _cur=savRun;
    return res;
}


