#include "../inputGenerator.h"

#include <iostream>
#include <fstream>

#include "../tools.h"
#include "../readInput.h"
#include "../multiIndex.h"
#include "../autoconverge.h"
#include "../problemTree.h"
#include "../mpiWrapper.h"
#include "printOutput.h"
#include "runDir.h"

std::shared_ptr<InputGenerator> InputGenerator::factory(int argc, char **argv){
    std::shared_ptr<InputGenerator> res,res1,res2,res3,res4;
    std::vector<std::shared_ptr<InputGenerator>> ress;
    ReadInput inp("",argc,argv);
    if(argc<2)exit(0);
    ress.push_back(std::shared_ptr<InputGenerator>(new Autoconverge(inp)));
    ress.push_back(std::shared_ptr<InputGenerator>(new InputRanges(argc,argv)));
    ress.push_back(std::shared_ptr<InputGenerator>(new InputSingle(argc,argv)));
    for(auto r: ress)if(r->nInps()>0)return r;
    ABORT("no input specified on command line");
    return res;
}

InputSingle::InputSingle(int argc, char **argv){
    // take input file "as is", last resort if all others do not generate
    _nxt=ReadInput::getInputFile("",argc,argv);
}

std::vector<std::string> InputGenerator::readInpFile(std::string File){
    if(not folder::exists(File))ABORT("cannot find input file "+File);

    std::ifstream inpc(File.c_str());
    std::vector<std::string> res(1);
    while(inpc.good()){
        std::getline(inpc,res.back());
        res.push_back("");
    }
    inpc.close();
    res.pop_back();
    return res;

}

std::string InputGenerator::writeInpFile(std::string Root,const std::vector<std::string> &InputLines){
    std::string res;
    if(MPIwrapper::isMaster(MPIwrapper::worldCommunicator())){
        if(not folder::exists(Root))folder::create(Root);

        res=_reuse->match(InputLines);

        // update run-directory (if changed)
        _runDir.updateRoot(Root);

        // we may need more digits for numbering all runs
        _runDir.updateDigits(std::max(4,int(std::log10(_runDir.nextUnused().cur()+nInps())+2)));

        if(res==""){
            std::string outDir=_runDir.nextUnused().dir();
            if(not folder::create(outDir))ABORT("could not create folder "+outDir);
            std::ofstream inpc((outDir+ReadInput::inputCopy).c_str());
            for(auto line: InputLines)inpc<<line<<std::endl;
            inpc.close();
            res=outDir+ReadInput::inputCopy;
        }
    }
    MPIwrapper::Bcast(res,MPIwrapper::master());
    return res;
}

InputRanges::InputRanges(int argc, char *argv[]):_current(-1)
{
    ReadInput inp(argc,argv);
    _reuse.reset(new ReuseRun(inp));

    std::string file(ReadInput::getInputFile("",argc,argv));
    std::ifstream stream(file.c_str(),std::ios::in);
    if(not stream.is_open())ABORT("could not open input file '"+file+"'");

    _root=file;
    for(auto e: {ReadInput::inputExtension,"/"+ReadInput::inputCopy})_root=tools::stringStrip(_root,e);

    std::vector<std::string> allLines;
    std::string line;
    while(std::getline(stream,line))allLines.push_back(line);
    stream.close();

    // find RANGE{...} macros in the lines (if any)
    std::vector<std::string> ranges;
    for(auto line: allLines){
        line=line.substr(0,line.find_first_of(ReadInput::comments));
        size_t p(0);
        while(std::string::npos!=(p=line.find("RANGE{"))){
            ranges.push_back(line.substr(p,line.find("}")-p+1));
            line=line.substr(p+ranges.back().length());
        }
    }
    if(ranges.size()==0)return; // no ranges specified

    if(file!=_root+ReadInput::inputExtension)
        ABORT("only files with extension "+ReadInput::inputExtension+" may contain RANGE{...} macros");

    std::vector<int> dims;
    std::vector<std::vector<std::string>> pars;
    for(size_t k=0;k<ranges.size();k++){
        auto r=tools::stringInBetween(ranges[k],"RANGE{","}",true);
        if(r.find(":")!=std::string::npos){
            auto grid=tools::rangeToGrid(r,0,{"",":",""});
            r="";
            for(auto g: grid)r+=tools::str(g)+",";
            r.pop_back();
        }
        pars.push_back(tools::splitString(r,','));
        dims.push_back(pars.back().size());
    }
    MultiIndex idx(dims);
    std::vector<int> next;
    while(idx.next(next)){

        std::vector<std::string>inpLines=allLines;
        for(size_t k=0;k<ranges.size();k++){
            for(auto &l: inpLines)l=tools::substringReplaceAll(l,ranges[k],pars[k][next[k]]);
            inpLines.insert(inpLines.begin(),std::string("#"+ranges[k]+" "+pars[k][next[k]]));
        }
        _inpLines.push_back(inpLines);
    }
}

std::string InputRanges::nextInput(){
    return (++_current)<int(_inpLines.size())?writeInpFile(_root,_inpLines[_current]):"";
}
int InputRanges::nInps() const{
    return _inpLines.size();
}

