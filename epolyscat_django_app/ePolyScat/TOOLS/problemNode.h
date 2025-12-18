#ifndef PROBLEMNODE_H
#define PROBLEMNODE_H

#include <vector>
#include <string>
#include <memory>
#include <fstream>

#include "readInputList.h"
#include "tools.h"

class ReadInputList;

class ProblemObservable;
class ProblemName;

/// abstract base class for nodes of a problem tree
class ProblemNode
{
    static std::string _sep;
    std::string _kind;
protected:
    /// split Definition, removing the derived class's kind(), return empty if not kind()
    std::vector<std::string> definitionSplit(const std::string Definition);

    /// join kind() and Parts into standard definition of node
    std::string definitionJoin(const std::vector<std::string> & Parts) const;

    std::vector<std::string> _inputLines;
    // return all lines of a given Category
    void allLinesInCategory(std::string InpcFile, std::string Category, std::vector<std::string> & Lines );
    // dump _inputLines into Stream
    void writeInputLines(std::ofstream & Stream) const{for(auto l: _inputLines)Stream<<l<<std::endl;}

    std::vector<std::string> inputItem(ReadInputList & Inp, std::string Category, std::string Name);
    void appendSingleInputLine(ReadInputList &Inp, std::string Category, std::string Name, std::vector<std::string> &InputLines);
public:
    // create derived classes from Definition string
    static std::shared_ptr<ProblemNode> factory(std::string Definition);
    virtual ~ProblemNode()=default;
    ProblemNode(std::string Kind):_kind(Kind){}

    virtual std::string definition() const =0;
    virtual void addToInput(std::ofstream &Stream) const=0;
    bool isDefinition(std::string String ) const {return String.find(_sep)!=std::string::npos;}

    virtual bool operator==(const ProblemNode & Other) const {return definition()==Other.definition();}
    virtual bool operator<(const ProblemNode &Other) const  {return definition()<Other.definition();}

    virtual void addToHeader(std::ofstream & Stream) const{Stream<<"";/* suppress unused variable warning */ };

    const std::string &kind() const {return _kind;}

    // human-readable info
    virtual std::string str() const {return definition();}

};
#endif // PROBLEMNODE_H
