#include "problemName.h"
#include "abort.h"
#include "str.h"

#include "problemTemplates.h"


ProblemName::ProblemName(ReadInputList &Inp):ProblemName()
{
    auto res=inputItem(Inp,"Problem","name");
    if(res.size()==1){
        if(_name=="")_name=res[0];
        if(res[0]!=_name)ABORT("conflicting Problem names on file"+res[0]+" and command line input="+_name);
    }
    if(_name=="")_name="ProblemAutocreate";

    ReadInput::main.read("Problem","description",_description,"","one-line description of the problem",1,"describeProblem");

     allLinesInCategory(Inp.outputTopDir()+"inpc",ReadInput::macro,_inputLines);
}
void ProblemName::addToHeader(std::ofstream &Stream) const{
    Stream<<std::endl;
    Stream<<"#     "+std::string(5,'+')+"  "+_name+" --- "+_description+"  "+std::string(5,'+')<<std::endl;
}

void ProblemName::addToInput(std::ofstream &Stream) const{
    Stream<<std::endl;
    Stream<<"Title:"<<std::endl;
    Stream<<_name<<" --- "<<_description<<"\n\n"<<std::endl;


    Stream<<std::endl;
    Stream<<"#Problem"+ReadInput::categoryTerminator+" name="+_name<<std::endl;
    Stream<<"#Problem"+ReadInput::categoryTerminator+" description="+_description<<std::endl;

    // add macros (if any)
    Stream<<std::endl;
    for(auto l: _inputLines)
        Stream<<l<<std::endl;
}

