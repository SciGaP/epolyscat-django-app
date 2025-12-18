#include "problemState.h"
#include "str.h"

ProblemState::ProblemState(ReadInputList &Inp):ProblemState()
{
    appendSingleInputLine(Inp,"Eigen","select",_inputLines);
    appendSingleInputLine(Inp,"Operator","initial",_inputLines);
    appendSingleInputLine(Inp,"Initial","kind",_inputLines);
    appendSingleInputLine(Inp,"Initial","state",_inputLines);
    _inputLines.size()>0?
                _inputLines.insert(_inputLines.begin(),"# --- state(s) used in calculation")
              : _inputLines.insert(_inputLines.begin(),"# --- default state(s) used");

}


void ProblemState::addToInput(std::ofstream &Stream) const{
    Stream<<std::endl;
    writeInputLines(Stream);
}
