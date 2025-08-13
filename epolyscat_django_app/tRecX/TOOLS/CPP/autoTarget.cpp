#include "autoTarget.h"
#include "str.h"

std::vector<std::vector<double>> AutoTarget::cols(std::string Dir){
    _fileKind=_def.substr(0,_def.find("["));
    if(_fileKind=="default"){
        for(auto d: _defaults){
            if(folder::exists(Dir+d.substr(0,d.find("["))))_fileKind=d.substr(0,d.find_first_of("{["));
        }
    }

    std::string def;
    for(auto d: _defaults)
        if(_fileKind==d.substr(0,d.find("[")))def=d;
    if(def=="")def="[all][all]";

    if(not folder::exists(Dir+file()))return {};

    _rows=tools::stringInBetween(_def.find("][")!=std::string::npos?_def:def,"[","][",true);

    // get columns (or default)
    std::string scol=_def.rfind("[")!=std::string::npos?_def:def;
    scol=scol.substr(scol.rfind("[")+1);
    scol=scol.substr(0,scol.rfind("]"));
    std::vector<int> cc=tools::string_to_intVec("{"+scol+"}");
    _cols.clear();
    for(auto c: cc)_cols.push_back(c);



    // read cols
    AsciiFile f(Dir+file());
    auto pars=rowPars();

    std::vector<std::vector<double>> ref(1);
    if(pars.size()==3)f.readCols(ref,{(unsigned int) pars.back()});

    std::vector<std::vector<double>> res;
    f.readCols(res,_cols);

    if(res[0].size()==0)return {};

    // select rows
    for(int k=res[0].size()-1;k>=0;k--){
        if(not inTarget(k,ref[0],pars)){
            for(size_t l=0;l<res.size();l++){
                res[l].erase(res[l].begin()+k);
            }
        }
    }
    return res;
}

std::string AutoTarget::str() const{
    auto pars=rowPars();
    std::string s("File: "+file());
    if(pars.size()==0)s+=",  Rows: all";
    else if(pars.size()==2)s+=Str(",  Rows: ","")+int(pars[0])+"-"+int(pars[1]);
    else if(pars.size()==3)s+=Str(" where col ","")+int(pars[2])+" in ["+pars[0]+","+pars[1]+"] ";
    s+=",  Cols:"+tools::str(_cols);
    return s;
}


std::vector<double> AutoTarget::rowPars() const {
    if(_rows=="all")return {};
    std::vector<double> res;
    if(_rows.find(":")!=std::string::npos){
        res.push_back(tools::string_to_double(tools::splitString(_rows,':')[0]));
        res.push_back(tools::string_to_double(tools::splitString(_rows,':')[1]));
    }
    if(_rows.find(",")!=std::string::npos){
        auto parts=tools::splitString(_rows,',');
        res.push_back(tools::string_to_double(parts[0]));
        res.push_back(tools::string_to_double(parts[1]));
        res.push_back(parts.size()==3?tools::string_to_double(parts[2]):0);
    }

    if(res.size()==0)
        ABORT("need rows as integer range [3:7] or value in reference column [1.5,3.7,2] for 2 in [1.5,3.7], got: {"+_rows+"}");
    return res;
}

bool AutoTarget::inTarget(int K, std::vector<double> RefCol, std::vector<double> Pars){
    if(Pars.size()==0)
        return true;
    if(RefCol.size()>0){
        return Pars[0]<=RefCol[K] and RefCol[K]<= Pars[1];
    }
    if(Pars.size()==2){
        return Pars[0]<=K and K<=Pars[1];
    }
    DEVABORT("inTarget case not covered");
    return false;
}



