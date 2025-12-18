#include "problemParameters.h"

#include <fstream>
#include "tools.h"
#include "autoconverge.h"

#include "printOutput.h"

ProblemParameters::ProblemParameters(ReadInputList &Inp):ProblemParameters()
{
    std::string fil=Inp.outputTopDir()+"conv",lin;
    std::ifstream f(fil.c_str());
    while(f.good()){
        std::getline(f,lin);
        if(lin.find("Parameter:")==0){
            std::vector<std::string> parts=tools::splitString(lin.substr(lin.find(":")+1),',');
            _cats.push_back(tools::cropString(parts[0]));
            _nams.push_back(tools::cropString(parts[1]));
            _lins.push_back(tools::cropString(parts[2]));
        }
    }
}

void ProblemParameters::addToInput(std::ofstream &Stream) const{
    Stream<<"";
}

void ProblemParameters::addToHeader(std::ofstream &Stream) const{
    Stream<<std::endl;

    Stream<<"#Autoconverge: category, name, atLine, lowerLimit, upperLimit, step"<<std::endl;
    Stream<<"#"+std::string(6,'-')+" uncomment the following lines to re-calculate error estimate ---"<<std::endl;
    Stream<<"#"<<std::endl;
    //    Stream<<"#"<<_observable->targetAndMetric()<<std::endl;
    Stream<<"#Autoconverge: category, name, atLine, lowerLimit, upperLimit, step"<<std::endl;
    for(size_t k=0; k<_cats.size();k++)
        Stream<<"# "<<_cats[k]<<", "<<_nams[k]<<", "<<_lins[k]<<std::endl;
    Stream<<"#"<<std::endl;
    Stream<<"#"+std::string(6,'-')+" end uncomment "+std::string(40,'-')<<std::endl;
}

std::string ProblemParameters::str() const{
    std::string s;
    for(size_t k=0;k<_cats.size();k++){
        s+=(_cats[k]!=Autoconverge::macro?
                    _cats[k]+":"+_nams[k]+"["+_lins[k]+"], "
                  :_nams[k]+", ");
    }
    if(s=="")return "Parameters: -- NONE varied --";
    s.pop_back();
    s.pop_back();
    return "Convergence: wrt. "+s;

}
