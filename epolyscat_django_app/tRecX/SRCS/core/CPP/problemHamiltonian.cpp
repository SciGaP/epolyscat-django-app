#include "../problemHamiltonian.h"

ProblemHamiltonian::ProblemHamiltonian(ReadInputList &Inp):ProblemNode("Hamiltonian")
{
    std::vector<std::string> inps;
    inps=inputItem(Inp,"Operator","hamiltonian");
    if(inps.size()==1)_hamiltonian=inps.front();
    inps=inputItem(Inp,"Operator","interaction");
    if(inps.size()==1)_interaction=inps.front();
}


void ProblemHamiltonian::addToInput(std::ofstream &Stream) const{
    Stream<<std::endl;
    Stream<<"Operator"+ReadInput::categoryTerminator+" hamiltonian="+_hamiltonian<<std::endl;
    if(_interaction!="")
        Stream<<"Operator"+ReadInput::categoryTerminator+" interaction="+_interaction<<std::endl;
}
