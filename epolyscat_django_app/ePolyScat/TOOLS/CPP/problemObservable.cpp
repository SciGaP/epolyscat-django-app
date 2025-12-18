#include "problemObservable.h"

#include "tools.h"
#include "str.h"

ProblemObservable::ProblemObservable(ReadInputList &Inp):ProblemObservable()
{
    std::string fil=Inp.outputTopDir()+"conv",lin;
    std::ifstream f(fil.c_str());
    while(f.good()){
        std::getline(f,lin);
        if(lin.find("Target:")==0){
            lin=lin.substr(lin.find(':')+1);
            _kind=tools::cropString(lin.substr(0,lin.find(':')));
            std::vector<std::string> parts=tools::splitString(lin.substr(lin.find(':')+1),',');
            for(auto p: parts){
                if(p.find("File:")!=std::string::npos)_file=tools::cropString(p.substr(p.find("File:")+5));
                if(p.find("Rows:")!=std::string::npos)_rows=tools::cropString(p.substr(p.find("Rows:")+5));
                if(p.find("Cols:")!=std::string::npos)_cols=tools::cropString(p.substr(p.find("Cols:")+5));
            }
        }
        if(lin.find("Metric:")==0)
            _metric=tools::cropString(lin.substr(lin.find(":")+1));
    }
}

// explanatory names for internal codes
static std::map<std::string,std::string> aka={
    {"eig","Eigenvalues"},
    {"RelativeDistanceWeight","L2-norm of element-wise relative accuracies of observable vector(s)"}
};

void ProblemObservable::addToHeader(std::ofstream &Stream) const{
    Stream<<"#"<<std::endl;
    Stream<<"#     observable: "+aka[_file]+", Rows: "+_rows+", Cols: "+_cols<<std::endl;
    Stream<<"#    errorMetric: "+aka[tools::splitString(_metric,',')[0]]+" "+_metric.substr(_metric.find(',')+1)<<std::endl;

}
std::string ProblemObservable::targetAndMetric() const{
    return "Autoconverge: target="+_file+"["+_rows+"]["+_cols+"]";
}


std::string ProblemObservable::str() const {
    std::string s;
    s=s+aka[_file]+", Rows: "+tools::str(_rows)+", Cols: "+tools::str(_cols);
    return "Observable: "+s;
}
