#ifndef PROBLEMNAME_H
#define PROBLEMNAME_H

#include "problemNode.h"
#include "problemTemplates.h"


class ProblemName : public ProblemNode
{
    std::string _name;
    std::string _description;

public:
    ProblemName(ReadInputList & Inp);

    ProblemName():ProblemNode("Name"){}
    std::string definition() const {
        if(isDefinition(_name))ABORT("bad: "+_name);
        std::vector<std::string> def({_name,_description});
        def.insert(def.end(),_inputLines.begin(),_inputLines.end());
        return definitionJoin(def);}
    ProblemName(std::string Definition)
        :ProblemName(){
        auto parts=definitionSplit(Definition);
        if(parts.size()>1){
            _name=parts[0];
            _description=parts[1];
            _inputLines={parts.begin()+2,parts.end()};
        }

    }

    void addToInput(std::ofstream &Stream) const;
    void addToHeader(std::ofstream &Stream) const;

    bool operator<(const ProblemName &Other) const {return _name.compare(Other._name);}
    std::string str() const {return _name;}
};

#endif // PROBLEMNAME_H
