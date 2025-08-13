#ifndef PROBLEMDATABASE_H
#define PROBLEMDATABASE_H

#include "problemNode.h"

class ProblemDatabase : public ProblemNode
{
    std::string _path;
public:
    ProblemDatabase():ProblemNode("Database"){}
    ProblemDatabase(std::string Definition):ProblemDatabase(){
        if(isDefinition(Definition)){
            auto parts=definitionSplit(Definition);
            if(parts.size()==1)_path=parts[0];
        }
        else
            _path=Definition;
    }

    std::string definition() const {return definitionJoin({_path});}
    void addToInput(std::ofstream &Stream) const {Stream<<"";};
    void addToHeader(std::ofstream &Stream) const {Stream<<"# --- from Database: "+_path+"\n";};
    std::string str() const {return "Database: "+_path;}
};

#endif // PROBLEMDATABASE_H
