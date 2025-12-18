#ifndef PROBLEMCOORDINATES_H
#define PROBLEMCOORDINATES_H

#include "problemNode.h"
#include "problemTemplates.h"

class ProblemCoordinates : public ProblemNode
{
    std::string _coors;
public:
    ProblemCoordinates(ReadInputList& Inp);

    ProblemCoordinates():ProblemNode("Coordinates"){}
    std::string definition() const {return definitionJoin({_coors});}
    ProblemCoordinates(std::string Definition):ProblemCoordinates(){
        auto parts=definitionSplit(Definition);
        if(parts.size()==1)_coors=parts[0];
    }

    void addToInput(std::ofstream &Stream) const;

    bool operator==(const ProblemNode & Other) const {
        auto o=dynamic_cast<const ProblemCoordinates*>(&Other);
        return o and not ((*this)<(*o)) and not ((*o)<(*this));
    };
    bool operator<(const ProblemCoordinates &Other) const {return _coors.compare(Other._coors);}
    std::string str() const {return "Coordinates: "+_coors;}
};

#endif // PROBLEMCOORDINATES_H
