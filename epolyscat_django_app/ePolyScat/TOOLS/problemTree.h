#ifndef PROBLEMTREE_H
#define PROBLEMTREE_H

#include <vector>
#include <string>
#include <memory>
#include "tree.h"
#include "inputGenerator.h"

#include "problemNode.h"
#include "problemResult.h"

class ReadInput;
class ProblemTree: public Tree<ProblemTree>, public InputSingle
{

protected:
    typedef std::shared_ptr<ProblemNode> Node;
    Node _node;
    int defaultBranch() const{return 0;}; // default problem from this node
    void dump(std::ofstream &Stream) const;    // dump to some file (attach modification history)
    void extractProblem(); // guide through problem tree and create input file
    void addConverged(std::string Dir); // add converged runs from directory tree to problem database
    void presentChoices() const; // show current node value, return next branch
    int promptNode() const; // show current node value, return next branch
    virtual std::vector<Node> distillRun(std::string RunDir) {DEVABORT("not implemented");};
public:
    ProblemTree(){}
    ProblemTree(int argc, char* argv[]);
    ProblemTree(std::string Path); // retrieve from some file

    int nInps() const override {return int(bool(_node));}
    std::shared_ptr<ProblemNode> node() const {return _node;};
    void insert(std::vector<Node> Nodes);
    void write(std::string Path) const;    // dump to some file (attach modification history)
    std::shared_ptr<ProblemResult> retrieve(std::vector<std::string> PathToResults);
    std::string availableAtNode() const;
    std::string createInputFile(std::string Root) const;

    std::string strNode(int Level = Tree_defaultKind) const override{
        if(Level==Tree_defaultKind)
            return node()?tools::abbreviate(node()->definition(),40):"--empty--";
        else
            DEVABORT("Level value not implemented");
        return "";
    }


};

#endif // PROBLEMTREE_H
