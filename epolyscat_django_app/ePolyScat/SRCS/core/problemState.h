#ifndef PROBLEMSTATE_H
#define PROBLEMSTATE_H

#include "problemNode.h"
#include "stringTools.h"

class ProblemState : public ProblemNode
{
public:
    ProblemState():ProblemNode("State"){}
    std::string definition() const {return definitionJoin(_inputLines);}
    ProblemState(std::string Definition):ProblemState(){
        _inputLines=definitionSplit(Definition);
    }

    ProblemState(ReadInputList& Inp);

    void addToInput(std::ofstream &Stream) const;
    std::string str() const {
        return "State: "
                +(_inputLines.size()<2?" -- default states --"
                                     :tools::str(std::vector<std::string>({_inputLines.begin()+1,_inputLines.end()})));}

};

#endif // PROBLEMSTATE_H
