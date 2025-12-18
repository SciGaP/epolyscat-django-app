#include "problemTree.h"

#include "readInput.h"
#include "printOutput.h"

#include "problemTemplates.h"

#include "problemName.h"
#include "problemObservable.h"
#include "problemResult.h"
#include "problemParameters.h"
#include "problemDatabase.h"
#include "runDir.h"


typedef std::shared_ptr<ProblemNode> Node;

ProblemTree::ProblemTree(int argc, char* argv[]){
    if(argc<2)return;
}

void ProblemTree::insert(std::vector<Node> Nodes){
    if(Nodes.size()==0)return; // noting to be inserted.

    // find matching branch
    ProblemTree* cTree(0);
    for(auto c: children()){
        if(c->node()->definition()==Nodes.front()->definition())cTree=c;
    }
    if(not cTree){
        // no matching branch, create new
        cTree=new ProblemTree();
        childAdd(cTree);
        cTree->_node=Nodes.front();
    }
    cTree->insert({Nodes.begin()+1,Nodes.end()});
}

std::string ProblemTree::createInputFile(std::string Root) const {
    // create a new input file in Root

    // if Prob is not leaf, descend default path from this Prob
    if(not isLeaf())
        return child(defaultBranch())->createInputFile(Root);

    // get path to this leaf
    std::vector<Node> path;
    for(auto p=this;p!=0;p=p->parent()){
        path.insert(path.begin(),p->node());
    }

    std::string fil=Root+ReadInput::inputExtension;
    std::ofstream f(fil.c_str());

    // create informative header in inpc
    f<<"# --- machine-generated input: any edit may invalidate error estimate ---\n\n"<<std::endl;
    if(kind<ProblemDatabase>(path))kind<ProblemDatabase>(path)->addToHeader(f);
    if(kind<ProblemName>(path))kind<ProblemName>(path)->addToHeader(f);
    if(kind<ProblemObservable>(path))kind<ProblemObservable>(path)->addToHeader(f);
    if(kind<ProblemResult>(path))kind<ProblemResult>(path)->addToHeader(f);
    if(kind<ProblemParameters>(path))kind<ProblemParameters>(path)->addToHeader(f);

    // write proper intput lines to inpc, as appropriate
    for(auto p: path){
        p->addToInput(f);
    }

    f.close();
    return fil;
}

void ProblemTree::write(std::string Path) const{
    std::ofstream stream(Path.c_str());
    dump(stream);
    PrintOutput::message("Problem tree on "+Path);
}

void ProblemTree::dump(std::ofstream &Stream) const{
    Stream<<node()->definition()<<std::endl;
    Stream<<childSize()<<std::endl;
    for(auto c: children())c->dump(Stream);
}

void ProblemTree::presentChoices() const{

    std::vector<const ProblemTree*> path;
    std::vector<int> indent;
    for(const ProblemTree* n=this;n!=0;n=n->parent()){
        path.insert(path.begin(),n);
        size_t pCol=path[0]->node()->str().find(":");
        if(pCol>30)pCol=-2;
        indent.insert(indent.begin(),pCol);
    }
    size_t maxColon=indent.size()==0?0:*std::max_element(indent.begin(),indent.end());

    for(size_t k=0;k<path.size();k++){
        std::cout<<std::string(maxColon-indent[k],' ')<<path[k]->node()->str()<<std::endl;
    }
}

int ProblemTree::promptNode() const{

    std::cout<<"--- problem definition ---"<<std::endl;
    presentChoices();

    // present choices
    std::cout<<"Next choices:"<<std::endl;
    for(size_t k=0;k<childSize();k++){
        std::cout<<"  ["<<k<<"] "<<child(k)->node()->str()<<std::endl;
    }
    std::cout<<"Enter number/all/quit: "<<std::flush;
    std::string ans;
    std::getline(std::cin,ans);
    if(ans=="")return defaultBranch();
    else if(ans[0]=='a')return -1;
    else if(ans[0]=='q')exit(0);
    int res=tools::string_to_int(ans);
    if((res==0 and tools::cropString(ans)!="0") or res>=int(childSize())){
        ABORT("not a valid input: "+ans+", allowed: 0..."+tools::str(childSize()-1)+" or a for all, or empty for default");
    }
    return res;
}

void ProblemTree::extractProblem(){
    if(not isLeaf()){
        int ans(0);
        if(childSize()>1)ans=promptNode();
        if(ans==-1)DEVABORT("not implemented yet");
        child(ans)->extractProblem();
    }
    else {
        PrintOutput::title("created model input file: "+createInputFile("extracted"));
        PrintOutput::paragraph();
        presentChoices();
        exit(0);
    }
}

void ProblemTree::addConverged(std::string Dir){
    std::vector<Node> path;
    if(folder::exists(Dir))path=distillRun(Dir);
    if(kind<ProblemResult>(path) and
            kind<ProblemResult>(path)->hasConvergence())insert(path);

    RunDir dir(Dir);
    while(dir.dir()!=""){
        if(folder::exists(dir.dir())){
            path=distillRun(dir.dir());
            if(kind<ProblemResult>(path) and kind<ProblemResult>(path)->hasConvergence())
                insert(path);
        }
    }
}


