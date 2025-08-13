#include "problemResultTrecX.h"

#include "tools.h"
#include "str.h"

#include "problemParameters.h"

ProblemResultTrecX::ProblemResultTrecX(ReadInputList &Inp)
{
    std::string lin;
    std::ifstream f;

    _dir=Inp.outputTopDir();

    // get results for the parameter
    f.open((_dir+"conv").c_str());
    while(f.good()){
        std::getline(f,lin);
        if(lin.find("Parameter:")==0)_res.push_back(res(lin));
    }
    f.close();

    // extract the host name and parallel processes
    f.open((_dir+"outf").c_str());
    int cnt=0;
    std::string hostMarker("Master host =");
    _host="--unknown--";
    _nproc=1;
    while(f.good() and cnt<100){
        cnt++;
        std::getline(f,lin);
        if(lin.find(hostMarker)!=std::string::npos){
            _host=lin.substr(lin.find(hostMarker)+hostMarker.length()+1);
            _host=_host.substr(0,std::min(_host.find("***"),_host.find(" (")));
            _host=tools::cropString(_host);
            if(lin.find("processes"))
                _nproc=tools::string_to_int(tools::stringInBetween(lin," (","processes)",true));
            break;
        }
    }
    f.close();

    // get total time for reference run
    f.open((_dir+"timer").c_str());
    std::string timeAll(" run_trecx: all ");
    _cpu=-1.;
    while(f.good()){
        std::getline(f,lin);
        // Note: luckily, split string crops string, if split char is ' '
        if(lin.find(timeAll)!=std::string::npos){
            _cpu=tools::string_to_double(tools::splitString(lin,' ')[3]);
        }
        if(_cpu>0.)break;
    }
    f.close();
}
