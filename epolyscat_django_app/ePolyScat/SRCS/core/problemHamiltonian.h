#ifndef PROBLEMHAMILTONIAN_H
#define PROBLEMHAMILTONIAN_H

#include "problemNode.h"

class ReadInput;
class ProblemHamiltonian : public ProblemNode
{
    std::string _hamiltonian,_interaction;
public:
    ProblemHamiltonian():ProblemNode("Hamiltonian"){}
    ProblemHamiltonian(ReadInputList& Inp);

    std::string definition() const{return definitionJoin({_hamiltonian,_interaction});}
    ProblemHamiltonian(std::string Definition)
        :ProblemHamiltonian(){
        auto parts=definitionSplit(Definition);
        if(parts.size()==2){_hamiltonian=parts[0],_interaction=parts[1];}

    }
    void addToInput(std::ofstream &Stream) const;

    bool operator<(const ProblemHamiltonian &Other) const {return (_hamiltonian+_interaction).compare(Other._hamiltonian+Other._interaction);}
    std::string str() const {return "Hamiltonian: "+_hamiltonian+(_interaction!=""?" + "+_interaction:"");}
};

#endif // PROBLEMHAMILTONIAN_H
