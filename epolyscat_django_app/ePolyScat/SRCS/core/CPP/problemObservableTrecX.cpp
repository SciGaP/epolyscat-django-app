#include "problemObservableTrecX.h"

#include "readInputList.h"

ProblemObservableTrecX::ProblemObservableTrecX(ReadInputList &Inp)
    :ProblemObservable(Inp)
{
    allLinesInCategory(Inp.outputTopDir()+"/inpc","Spectrum",_inputLines);
    allLinesInCategory(Inp.outputTopDir()+"/inpc","Eigen",_inputLines);
}
