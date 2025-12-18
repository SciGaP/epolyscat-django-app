#include "problemResult.h"

#include "tools.h"
#include "str.h"

#include "problemParameters.h"
#include "problemObservable.h"
#include "outputTable.h"

ProblemResult::res::res(std::string Line){
    std::vector<std::string> parts=tools::splitString(Line.substr(Line.find(":")+1),',');
    std::string lin=tools::cropString(parts[2]);
    _ref=tools::string_to_double(parts[3]);

    _err=parts[4];
    _def={parts.begin(),parts.begin()+3};
    for(size_t k=5;k<parts.size();k++){
        _pars.push_back(tools::string_to_double(tools::splitString(parts[k],':')[0]));
        _errs.push_back(tools::string_to_double(tools::splitString(parts[k],':')[1]));
    }
    if(ref()==(int)_pars.size())
        DEVABORT(Sstr+"reference parameter "+_ref+" not in list: "+_pars)
}

std::string ProblemResult::res::line() const {
    std::vector<std::string> parts(_def);
    parts[0]="Parameter: "+parts[0];
    parts.push_back(tools::str(_ref));
    parts.push_back(_err);
    auto pE=parsErrs();
    parts.insert(parts.end(),pE.begin(),pE.end());
    return tools::joinString(parts,",");
}

bool ProblemResult::hasConvergence() const{
    return _res.size()>1;
}

double ProblemResult::errTotal() const{
    double errSq=0;
    for(size_t k=0;k<_res.size();k++)
        errSq+=std::pow(_res[k]._errs[_res[k].ref()],2);
    return sqrt(errSq);
}

void ProblemResult::addToHeader(std::ofstream &Stream) const{

    OutputTable tab;
    tab.addRow({"Paramter[Line]","Value","ErrorEstimate"});
    for(size_t k=0;k<_res.size();k++){
        int n=_res[k].ref();
        tab.newRow();
        tab.rowItem(_res[k].parnam());
        tab.rowItem(_res[k]._pars[n]);
        tab.rowItem(_res[k]._errs[n]);
    }

    Stream<<std::endl;
    Stream<<"\n"
            "# --- WARNING --- error(s) below are only valid for specified observable ------\n"
            "#\n"
            "# conclusions about any other observables are perilious\n"
            "# error estimates are based on variation of the indicated parameters\n"
            "# these are likely to be the most critical for the given problem\n"
            "#\n"
            "#       --- no other inputs were examined for errors ---\n"
            "#\n"
            "# !!! any edit of this model input file will change the actual error !!!\n"
            "#-----------------------------------------------------------------------------\n"
            "#"<<std::endl;
    Stream<<"# Overall error: "<<sqrt(errTotal())<<std::endl;
    Stream<<"#           CPU: "<<_cpu<<" on Host "<<_host<<(_nproc>1?"("+tools::str(_nproc)+" threads) ":" ")<<_dir<<std::endl;
    Stream<<"#"<<std::endl;

    auto lines=tab.format();
    for(auto l: lines)Stream<<"# "<<l<<std::endl;

    Stream<<std::endl;

}
