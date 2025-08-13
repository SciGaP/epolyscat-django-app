#ifndef PROBLEMOBSERVABLE_H
#define PROBLEMOBSERVABLE_H

#include "problemNode.h"

class ProblemObservable : public ProblemNode
{
    std::string _kind;
    std::string _file;
    std::string _rows;
    std::string _cols;
    std::string _metric;
public:
    ProblemObservable():ProblemNode("Observable"){}
    std::string definition() const {
        return definitionJoin({_kind,_file,_rows,_cols,_metric});}
    ProblemObservable(std::string Definition):ProblemObservable(){
        auto parts=definitionSplit(Definition);
        if(parts.size()==5){
            _kind=parts[0];
            _file=parts[1];
            _rows=parts[2];
            _cols=parts[3];
            _metric=parts[4];
        }
    }

    ProblemObservable(ReadInputList& Inp);



    void addToInput(std::ofstream &Stream) const{ writeInputLines(Stream);}
    void addToHeader(std::ofstream &Stream) const;
    std::string targetAndMetric() const;

    bool operator==(const ProblemNode & Other) const {
        auto o=dynamic_cast<const ProblemObservable*>(&Other);
        return o and not ((*this)<(*o)) and not ((*o)<(*this));
    };
    bool operator<(const ProblemObservable &Other) const {
        return _file.compare(Other._file);
        return _cols.compare(Other._cols);
        return _rows.compare(Other._rows);
    }
    std::string str() const;
};

#endif // PROBLEMOBSERVABLE_H
