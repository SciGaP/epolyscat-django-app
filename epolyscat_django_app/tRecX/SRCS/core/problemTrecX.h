#ifndef PROBLEMTRECX_H
#define PROBLEMTRECX_H


#include "problemTree.h"
#include "problemNode.h"

typedef std::shared_ptr<ProblemNode> Node;


class ProblemTrecX : public ProblemTree
{
    ProblemTrecX(std::ifstream &Stream);
    void loadDatabase(std::string Path);
    static std::shared_ptr<ProblemNode> nodeFactory(std::string Def);
    static bool readNode(std::ifstream &Stream,size_t & ChildSize, Node & Node);
    /// distill database entry from run-directory
    std::vector<Node> distillRun(std::string RunDir);
public:

    ProblemTrecX(int argc, char* argv[]);

};

#endif // PROBLEMTRECX_H
