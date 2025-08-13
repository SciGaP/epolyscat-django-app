#include "problemTrecX.h"

#include "abort.h"

#include "problemNode.h"
#include "problemName.h"
#include "problemHamiltonian.h"
#include "problemState.h"
#include "problemLaser.h"
#include "problemCoordinates.h"
#include "problemParameters.h"
#include "problemAxes.h"

#include "problemResultTrecX.h"
#include "problemObservableTrecX.h"
#include "problemDatabase.h"

typedef std::shared_ptr<ProblemNode> Node;

static std::string defaultDataBase="ProblemDataBase/test.db";

/// for -database=path_to_db_file creates database from file
/// for -addProblem=Run/Dir/1234 adds run to problem tree
ProblemTrecX::ProblemTrecX(int argc, char* argv[]):ProblemTree(argc,argv){
    std::string runDir,db=defaultDataBase;

    bool extract=false;
    for(int k=1;k<argc;k++){
        std::string argK(argv[k]);
        if(argK.find("-addProblem=")==0)runDir=tools::cropString(argK.substr(argK.find("=")+1));
        if(argK.find("-database=")==0)db=tools::cropString(argK.substr(argK.find("=")+1));
        if(argK=="-extractProblem")extract=true;
    }
    // no ProblemTree actions
    if(not extract and runDir=="")return;

    loadDatabase(db);

    if(extract)extractProblem();

    addConverged(runDir);


    write("ProblemDataBase/test.db");

    ABORT("created input file "+createInputFile("scrProb"));
}

std::vector<Node> ProblemTrecX::distillRun(std::string RunDir){
    ReadInputList inp(RunDir+"/"+ReadInput::inputList);

    // construct path - sequence of items determines tree hierarchy
    std::vector<Node> res;

    // problem definition
    res.push_back(Node(new ProblemName(inp)));
    res.push_back(Node(new ProblemHamiltonian(inp)));
    res.push_back(Node(new ProblemState(inp)));
    res.push_back(Node(new ProblemLaser(inp)));
    res.push_back(Node(new ProblemObservableTrecX(inp)));

    // discretization
    res.push_back(Node(new ProblemCoordinates(inp)));

    // parameters varied for the error estimate
    res.push_back(Node(new ProblemParameters(inp)));
    res.push_back(Node(new ProblemAxes(inp)));

    // results of the specific calculation
    res.push_back(Node(new ProblemResultTrecX(inp)));

    return res;
}

std::shared_ptr<ProblemNode> ProblemTrecX::nodeFactory(std::string Def){
    std::shared_ptr<ProblemNode> res;
    if((res=ProblemNode::factory(Def)))return res;
    if((res=makeNode<ProblemHamiltonian>(Def)))return res;
    if((res=makeNode<ProblemLaser>(Def)))return res;
    if((res=makeNode<ProblemCoordinates>(Def)))return res;
    if((res=makeNode<ProblemAxes>(Def)))return res;
    if((res=makeNode<ProblemState>(Def)))return res;
    if((res=makeNode<ProblemResult>(Def)))return res;
    return res;
}

// once functionality is good, switch to binary files
bool ProblemTrecX::readNode(std::ifstream &Stream,size_t & ChildSize, Node & Node){
    //    size_t len;
    //    Stream>>len;
    //    std::string def(len,' ');
    //    Stream>>def;
    std::string lin;
    std::getline(Stream,lin);
    Node=nodeFactory(lin);
    if(not Node)ABORT("unidentified: "+lin);
    std::getline(Stream,lin);
    ChildSize=tools::string_to_int(lin);
    return Stream.good();
}

void ProblemTrecX::loadDatabase(std::string DataBasePath){
    std::ifstream stream(DataBasePath.c_str());
    if(stream.good()){
        size_t cSize(0);
        if(not readNode(stream,cSize,_node))DEVABORT("not a database file "+DataBasePath);
        for(size_t k=0;k<cSize;k++)childAdd(new ProblemTrecX(stream));
    }
}


ProblemTrecX::ProblemTrecX(std::ifstream &Stream){
    size_t cSize(0);
    if(not readNode(Stream,cSize,_node))DEVABORT("database Stream corrupted");
    for(size_t k=0;k<cSize;k++)childAdd(new ProblemTrecX(Stream));
}



