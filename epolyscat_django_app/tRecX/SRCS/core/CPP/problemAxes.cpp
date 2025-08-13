#include "problemAxes.h"

#include "str.h"

ProblemAxes::ProblemAxes(ReadInputList &Inp):ProblemAxes()
{
    allLinesInCategory(Inp.outputTopDir()+"inpc","Axis",_inputLines);
}
void ProblemAxes::addToInput(std::ofstream &Stream) const{
    Stream<<std::endl;
    writeInputLines(Stream);
}

// we need this in class Axis, or so.
std::string axisBrief(std::vector<std::string> & Parts){
    for(auto &p: Parts)p=tools::cropString(p);
    std::string s;
    s+=Parts[0];
    if(Parts.size()>3)s+="["+Parts[2]+","+Parts[3]+"]";
    if(Parts.size()>1)s+=(s.back()==']'?"":"|")+Parts[1];
    if(Parts.size()>5)s+="@"+Parts[5];
    return s;
}

std::string ProblemAxes::str() const {
    std::string s;
    for(auto l: _inputLines){
        auto parts=tools::splitString(l,',');
        if(parts.size()==0 or parts[0].find("Axis:")==0)continue;
        s=s+"("+axisBrief(parts)+")x";
    }
    s.pop_back();
    return "Axes: "+s;
}
