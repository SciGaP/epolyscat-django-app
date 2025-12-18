#include "problemNode.h"

#include "abort.h"

#include "readInputList.h"
#include "tools.h"

#include "problemObservable.h"
#include "problemName.h"
#include "problemResult.h"
#include "problemParameters.h"
#include "problemDatabase.h"

std::string ProblemNode::_sep="-|-";

std::string ProblemNode::definitionJoin(const std::vector<std::string> &Parts) const{
    std::vector<std::string> parts(Parts);
    parts.insert(parts.begin(),kind());
    return tools::joinString(parts,_sep);
}
std::vector<std::string> ProblemNode::definitionSplit(const std::string Definition){
    auto res=tools::splitString(Definition,_sep);
    if(res.front()==kind())return {res.begin()+1,res.end()};
    else                   return{};
}

std::shared_ptr<ProblemNode> ProblemNode::factory(std::string Definition){
    std::shared_ptr<ProblemNode> res;
    if(bool(res=makeNode<ProblemDatabase>(Definition)))return res;
    if(bool(res=makeNode<ProblemName>(Definition)))return res;
    if(bool(res=makeNode<ProblemObservable>(Definition)))return res;
    if(bool(res=makeNode<ProblemResult>(Definition)))return res;
    if(bool(res=makeNode<ProblemParameters>(Definition)))return res;
    return res;
}

std::vector<std::string> ProblemNode::inputItem(ReadInputList &Inp, std::string Category, std::string Name){
    std::vector<std::string> res;
    do {
        res.push_back(Inp.readInput(Category,Name,res.size()+1));
    }
    while (res.back().find(ReadInput::notFound)==std::string::npos or res.back()=="");
    res.pop_back();
    return res;
}

void ProblemNode::appendSingleInputLine(ReadInputList &Inp, std::string Category, std::string Name, std::vector<std::string> &InputLines){
    auto its=inputItem(Inp,Category,Name);
    if(its.size()>0)InputLines.push_back(Category+ReadInput::categoryTerminator+" "+Name+"="+tools::str(its));
    if(its.size()>1)ABORT("cannot append single input line, found multiple inputs "+InputLines.back());
}


void ProblemNode::allLinesInCategory(std::string InpcFile, std::string Category, std::vector<std::string>& Lines){
    std::ifstream finp(InpcFile.c_str());
    bool inCat=false;
    while (finp.good()){
        std::string lin;
        std::getline(finp,lin);

        if(Category==ReadInput::macro and lin.find("#define ")==0)
            Lines.push_back(lin);

        // start accepting at Category bane
        if(lin.find(Category+ReadInput::categoryTerminator)==0)inCat=true;

        // stop accepting at new category or comment
        if(ReadInput::comments.find(lin[0])!=std::string::npos)inCat=false;
        if(ReadInput::isCategoryLine(lin) and lin.find(Category)!=0)inCat=false;
        if(inCat){
            Lines.push_back(lin);
            // stop accepting after blank line
            if(tools::cropString(lin)=="")inCat=false;
        }
    }
}

