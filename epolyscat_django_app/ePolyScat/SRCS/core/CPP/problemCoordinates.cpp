#include "problemCoordinates.h"

#include <algorithm>

#include "tools.h"
ProblemCoordinates::ProblemCoordinates(ReadInputList &Inp):ProblemCoordinates()
{
    auto coors=inputItem(Inp,"Axis","name");
    // remove duplicates
    std::sort(coors.begin(),coors.end());
    for(size_t k=coors.size();k>1;k--)
        if(coors[k-1]==coors[k])coors.erase(coors.begin()+k);
    _coors=tools::joinString(coors,".");
}
void ProblemCoordinates::addToInput(std::ofstream &Stream) const{
    writeInputLines(Stream);
}
