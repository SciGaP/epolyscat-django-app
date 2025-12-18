#include "autoconverge.h"
#include "autoTarget.h"
#include "targetError.h"

#include "readInput.h"
#include "printOutput.h"
#include "farm.h"
#include "asciiFile.h"
#include "metric.h"


const std::string Autoconverge::macro="_MACRO_";
std::string Autoconverge::estimate="EstimateError";

Autoconverge::Autoconverge(ReadInput &Inp, std::vector<Autoconverge*> Path)
{
    _reuse.reset(new ReuseRun(Inp,{"Autoconverge:"}));

    bool isAutoconverge=Inp.found("Autoconverge");

    // determine next line
    int line(Path.size()+1);
    std::string mode;
    Inp.read("Autoconverge","mode",mode,"ConvergeToTarget","Operation mode: "
                                                           "ConvergeToTarget...determine parameters for accuracy, "
                                                           "EstimateError...estimate error at given parameter");

    Inp.read("Autoconverge","accuracy",_accuracyTarget,"0","target relative accuracy");
    if(isAutoconverge and _accuracyTarget==0 and mode!=estimate){
        mode=estimate;
        if(Path.size()==0)PrintOutput::message("no accuracy target given, assuming mode="+mode);
    }
    if(_accuracyTarget!=0 and mode==estimate)
        PrintOutput::warning("Autoconverge:accuracy will be ignored with action="+mode);

    std::string target;
    Inp.read("Autoconverge","target",target,"default","target definition, e.g. spec_partial[all][2] all rows in second column in spec_partial");
    if(Path.size()==0)_target.reset(new AutoTarget(target,{"spec_partial[all][0,2]","spec_total[all][0,1]","eig[0:0][0]"}));

    double metric;
    std::string preCat(macro),atLine,strInit,strStep;

    _conNode.reset(new ConvergenceNode(Inp,Path.size()?Path[0]:this,line));

    Inp.read("Autoconverge","metric",metric,"1e-2","pedestal in relative distance for error definition",line);
    _metric.reset(new RelativeDistanceWeighted(metric));

    // do not actually create run input
    if(not Inp.found("Autoconverge"))return;

    // all levels share root and first input (this should be changed to "referenceInput()" and "rootDir()"
    if(Path.size()==0){
        _rootInp=Inp.root();
        _reference=tools::stringStrip(writeInpFile(_rootInp,Inp.inputLines()),ReadInput::inputCopy);
    }

    std::string nxtNam;
    Inp.read("Autoconverge","name",nxtNam,"",ReadInput::doNotAddToItemTable,line+1);
    if(nxtNam=="")return;

    Path.push_back(this);
    childAdd(new Autoconverge(Inp,Path));
    Path.pop_back();

}

int Autoconverge::nInps() const {
    int siz=1;
    for(const Autoconverge* n=root();n!=0;n=n->nodeNext()){
        if(node()->size()==0)return 0;
        siz+=n->node()->size();
    }
    return siz;
}

bool Autoconverge::runCompleted(std::string RunDir) const {
    return _target->cols(RunDir).size()>0;
}

std::string Autoconverge::nextInput() {
    std::string inp;
    size_t p=10;
    for(;p>0 and (""==(inp=nextPriority(p)));p--);
    if(inp=="")print();
    return inp;
}


std::string Autoconverge::nextPriority(size_t Priority){
    if(converged())return "";
    if(node()->pars().front()<node()->parMin())ABORT(Sstr+"cannot have lowerLimit="+node()->parMin()+"above initial parameter"+node()->pars().front());
    if(node()->pars().front()>node()->parMax())ABORT(Sstr+"cannot have upperLimit="+node()->parMax()+"below initial parameter"+node()->pars().front());

    // try update on lower level
    std::string next;
    for(auto c: children()){
        if(""!=(next=c->nextPriority(Priority)))return next;
    }

    if(node()->inps().size()==0){
        // this has maximal priority
        node()->inps().push_back(reference());
    }
    else {
        if(Priority<1)return "";
        if(not addNextParameter(Priority))return "";

        if(node()->pars().size()>1 and node()->pars().back()==node()->pars()[node()->pars().size()-2]){
            // no new parameter
            return "";
        }
        node()->inputRows()[node()->row()][node()->col()]=tools::str(node()->pars().back());

        // write new input, add new file name and return
        std::vector<std::string> inpLines;
        for(auto row: node()->inputRows()){
            inpLines.push_back(tools::joinString(row,(row.size()>0 and row[0]==std::string("#define"))?" ":", "));
        }
        node()->inps().push_back(tools::stringStrip(writeInpFile(referenceRoot(),inpLines),ReadInput::inputCopy));
    }
    return node()->inps().back()+ReadInput::inputCopy;
}

std::string Autoconverge::strNode(int Level) const{
    if(Level!=Tree_defaultKind)ABORT("only default");
    Str s(""," ");
    s=s+node()->nam()+" "+node()->cat()+node()->parMin()+node()->parMax();
    return std::move(s);
}


static std::string parAndErr(double Par, double Err){
    return tools::str(Par)+"("+tools::str(Err,2)+")";
}

int Autoconverge::parameterRow() const{
    for(size_t k=0;k<node()->inputRows().size();k++){
        std::vector<std::string> Row=node()->inputRows()[k];
        if(Row.size()){
            if(node()->cat()==macro and Row[0]=="#define" and Row[1].find(node()->nam())!=std::string::npos)return k;
            if(Row[0].find(node()->cat()+ReadInput::categoryTerminator)==0
                    and Row[0].substr(Row[0].find(ReadInput::categoryTerminator)).find(node()->nam())!=std::string::npos)return k;
        }
    }
    return node()->row()+1; // parameterRow after value indicates not found
}
void Autoconverge::print(){
    if(not MPIwrapper::isMaster(MPIwrapper::worldCommunicator()))return;

    std::vector<double> pars,delts,gams,alfs,errs;
    double parPred;

    node()->nextParameter(pars,errs,parPred); // need to call this to have target-info complete
    node()->parametersAndErrors(pars,errs);

    if(parent()==0){
        PrintOutput::end();
        PrintOutput::title("Error estimates and target parameter");
        PrintOutput::paragraph();
        PrintOutput::subTitle("    "+target()->str());
        PrintOutput::subTitle(Sstr+"   Runs:"+node()->inps());
        PrintOutput::paragraph();
        PrintOutput::newRow();
        PrintOutput::rowItem("Category:name");
        PrintOutput::rowItem("predicted function");
        PrintOutput::rowItem("predicted(error)");
        PrintOutput::rowItem("parameter(error)'s");
    }

    PrintOutput::newRow();
    PrintOutput::rowItem(node()->cat()+ReadInput::categoryTerminator+node()->nam());
    PrintOutput::rowItem(node()->errorFunction()->algebra());

    int k=std::find(pars.begin(),pars.end(),node()->pars()[0])-pars.begin();
    PrintOutput::rowItem(parAndErr(pars[k],errs[k]));

    std::string s="";
    for(size_t k=0;k<errs.size();k++)s+=parAndErr(pars[k],errs[k])+" ";
    s.pop_back();
    PrintOutput::rowItem(s);
    for(auto c: children())c->print();

    if(parent()==0){
        PrintOutput::paragraph();
        std::string f=writeConvergence();
        PrintOutput::subTitle("    Convergence data on: "+f);
        f=writeResults(); // result summary for machine read
        PrintOutput::subTitle(" Convergence results on: "+f);
    }
}


void Autoconverge::getDeltas(std::vector<double>& CPars, std::vector<double> &Delta, std::vector<double> &Gamma, std::vector<double> &Alpha) const
{
    std::vector<double> IPars(node()->pars());
    std::vector<std::string> Inps(node()->inps());
    IPars.resize(Inps.size());
    tools::sortByKey(IPars,Inps);


    // determine convergence factors between subsequent and to larges parameter

    // check for completed runs:
    std::vector<size_t> comp;
    CPars.clear();
    for(size_t k=0;k<Inps.size();k++){
        if(target()->cols(Inps[k]).size()>0){
            comp.push_back(k);
            CPars.push_back(IPars[k]);
        }
    }

    Delta.clear();
    Gamma.clear();
    Alpha.assign(1,1e-1);
    for(size_t k=1;k<comp.size();k++){
        Delta.push_back(distance(Inps[comp.back()],Inps[comp[k-1]]));
        Gamma.push_back(distance(Inps[comp[k]],Inps[comp[k-1]]));
        if(Gamma.size()>1)Alpha.push_back(Gamma.back()/Gamma[Gamma.size()-2]);
    }

}

std::string Autoconverge::writeConvergence(){
    std::vector<double>cpars,delta;
    node()->getDeltas(cpars,delta);

    AsciiFile f(tools::stringStrip(reference(),"inpc")+"convData");
    std::vector<std::string> head;
    std::string inps;
    for(size_t k=0;k<cpars.size();k++)inps+="  "+node()->inps()[k];
    head.push_back(inps);
    head.push_back(node()->nam()+"         d_best");
    f.writeComments(head);
    f.writeCols({{cpars.begin(),cpars.end()-1},delta});

    return f.name();
}

bool Autoconverge::addNextParameter(size_t Priority){

    if(Priority >1 and node()->pars().size()==1){
        node()->pars().push_back(node()->pars().front()-node()->parStep());
    }
    else if(Priority>2 and node()->pars().size()==2){
        node()->pars().push_back(node()->pars().front()+node()->parStep());
    }
    else if(Priority>5 and *std::max_element(node()->pars().begin(),node()->pars().end())<node()->parMax()){
        std::vector<double> pars,errs;
        double dum;
        double nxt=node()->nextParameter(pars,errs,dum);
        node()->pars().push_back(int(nxt+0.9));
        if(std::count(node()->pars().begin(),node()->pars().end(),node()->pars().back())>1)return false;;
    }
    else
        return false;

    return true;
}


double Autoconverge::distance(std::string DirA, std::string DirB) const {

    std::vector<std::vector<double>> colsA(target()->cols(DirA));
    std::vector<std::vector<double>> colsB(target()->cols(DirB));

    double dsqu=0;
    for(size_t c=0;c<colsA.size();c++){
        dsqu+=std::pow(metric()->distance(colsA[c],colsB[c]),2);
    }
    return sqrt(dsqu/colsA[0].size());
}

std::string Autoconverge::writeResults() const {
    std::ofstream stream;
    std::string f=reference()+"conv";
    stream.open(f.c_str());
    root()->writeResults(stream);
    return f;
}

void Autoconverge::writeResults(std::ofstream & ResultStream) const {
    if(parent()==0){
        ResultStream<<"Metric: "+metric()->definition()<<std::endl;
        ResultStream<<"Target: "+target()->definition()<<std::endl;
        ResultStream<<"# ---- Category, Name, Line, RefPar, ErrorFunction, Parameter:Error pairs"<<std::endl;
    }
    std::vector<double>pars,delts,gams,alfs,errs;
    getDeltas(pars,delts,gams,alfs);
    errs=node()->errorFunction()->estimates(pars,delts);
    ResultStream<<"Parameter: "+node()->cat()+", "+node()->nam()+", "+tools::str(node()->row()-parameterRow())+", "+tools::str(node()->pars()[0])+", "
            +node()->errorFunction()->algebra();
    for(size_t k=0;k<pars.size();k++)ResultStream<<", "+tools::str(pars[k],4)+":"+tools::str(errs[k],2);
    ResultStream<<std::endl;

    for(auto c: children())c->writeResults(ResultStream);

    if(parent()==0)ResultStream.close();
}


















