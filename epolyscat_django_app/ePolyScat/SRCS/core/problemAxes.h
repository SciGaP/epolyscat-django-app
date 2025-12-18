#ifndef PROBLEMAXES_H
#define PROBLEMAXES_H

#include "problemNode.h"

class ReadInputList;
class ProblemAxes : public ProblemNode
{
public:
    ProblemAxes(ReadInputList& Inp);

    ProblemAxes():ProblemNode("Axes"){}
    std::string definition() const {
        return definitionJoin(_inputLines);
    }
    ProblemAxes(std::string Definition):ProblemAxes(){
        _inputLines=definitionSplit(Definition);
    }

    void addToInput(std::ofstream &Stream) const;
    bool operator<(const ProblemNode &Other) const {return ProblemNode::operator<(Other);}
    std::string str() const;

};

#endif // PROBLEMAXES_H
