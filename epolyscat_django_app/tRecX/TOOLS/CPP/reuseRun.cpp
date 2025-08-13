#include "reuseRun.h"

#include "readInput.h"
#include "tools.h"
#include "str.h"

ReuseRun::ReuseRun(ReadInput &Inp,std::vector<std::string> Disregard):_disregard(Disregard){
    _disregard.push_back({"Reuse:"}); // never consider the Reuse category when comparing input files

    Inp.read("Reuse","runs",_reuse,"",
             "runs to re-use if identical input - blank-separated runDirs or [all] or format [1,4,6-9]",1,"reuse");
    if (_reuse=="[all]")
        _reuse=Inp.root()+"/[0-9999]";
    else if(_reuse.find("[")==0)
        _reuse=Inp.root()+"/"+_reuse.substr(_reuse.find("["));
}

static std::vector<int> string_to_intVec(std::string Range){
    Range=tools::substringReplaceAll(Range,"-",":");
    std::vector<std::string> parts=tools::splitString(Range,',');
    std::vector<int> res,tmp;
    for(auto p: parts){
        tmp=tools::string_to_intVec(p.find(":")==std::string::npos?"{"+p+"}":p);
        res.insert(res.end(),tmp.begin(),tmp.end());
    }
    return res;
}

static bool exception(std::string Line, std::vector<std::string> Except, std::string &Category){

    // update current category
    if(ReadInput::isCategoryLine(Line))
        Category=Line.substr(0,Line.find(ReadInput::categoryTerminator)+1);
    if(Line.find("#define ")==0)
        Category="#define ";
    else if(Line.find_first_of(ReadInput::comments)==0 or Line=="")
        Category="";

    // not in an input category
    if(Category=="")return true;

    // check for specific exceptions
    for(auto x: Except)
        if(Line.find(x)==0 or Category.find(x)==0){
            return true;
        }

    return false;
}

static std::string equalInps(std::string File, const std::vector<std::string> & Lines,std::vector<std::string> Except){

    if(not folder::exists(File))
        return "";

    std::ifstream inpf(File.c_str());
    auto pLin=Lines.begin();
    std::string cLin;
    std::string catFil="";
    std::string catLin="";
    while(inpf.good() and pLin<Lines.end()){
        // advance both sets of lines until not exception
        while (pLin<Lines.end() and exception(*pLin,Except,catLin))pLin++;
        do {std::getline(inpf,cLin);} while (inpf.good() and exception(cLin,Except,catFil));

        if(pLin<Lines.end() and inpf.good() and *pLin!=cLin )return "";
        pLin++;
    }

    while(std::getline(inpf,cLin).good())
        if(not exception(cLin,Except,catLin))return "";

    while(pLin<Lines.end()){
        if(not exception(*pLin,Except,catFil))return "";
        pLin++;
    }
    return File;
}

std::string ReuseRun::match(const std::vector<std::string> &InputLines) const {
    std::string res;
    if(_reuse.find("[")!=std::string::npos){
        std::vector<int> runs=string_to_intVec(tools::stringInBetween(_reuse,"[","]"));
        std::string root=_reuse.substr(0,_reuse.find("["));

        for(auto r: runs)
            if(""!=(res=equalInps(root+tools::str(r,4,'0')+"/"+ReadInput::inputCopy,InputLines,_disregard))){
                return res;
            }
    }
    else{
        std::vector<std::string> fils=tools::splitString(_reuse,',');
        for(auto f: fils)
            if(""!=(res=equalInps(f,InputLines,_disregard)))
                return res;
    }
    return "";
}
