#ifndef PROBLEMPARAMETERS_H
#define PROBLEMPARAMETERS_H

#include "problemNode.h"

class ReadInputList;
class Algebra;

/// list of parameters that where varied in the problem
class ProblemParameters : public ProblemNode
{
    std::vector<std::string> _cats,_nams,_lins;
public:
    ProblemParameters(ReadInputList & Inp);


    ProblemParameters():ProblemNode("Parameters"){}
    std::string definition() const {
        std::vector<std::string> parts(_cats);
        parts.insert(parts.end(),_nams.begin(),_nams.end());
        parts.insert(parts.end(),_lins.begin(),_lins.end());
        return definitionJoin(parts);}
    ProblemParameters(std::string Definition)
        :ProblemParameters(){
        auto parts=definitionSplit(Definition);
        size_t q3=parts.size()/3;
        _cats={parts.begin()+0*q3,parts.begin()+1*q3};
        _nams={parts.begin()+1*q3,parts.begin()+2*q3};
        _lins={parts.begin()+2*q3,parts.begin()+3*q3};
    }


    std::vector<std::string> &cats(){return _cats;}
    std::vector<std::string> &nams(){return _nams;}
    std::vector<std::string> &lins(){return _lins;}

    void addToHeader(std::ofstream &Stream) const;
    void addToInput(std::ofstream &Stream) const;

    bool operator==(const ProblemNode & Other) const {
        auto o=dynamic_cast<const ProblemParameters*>(&Other);
        return o and not ((*this)<(*o)) and not ((*o)<(*this));
    };
    bool operator<(const ProblemParameters &Other) const {return _cats[0].compare(Other._cats[0]);}
    std::string str() const;

};

#endif // PROBLEMPARAMETERS_H
