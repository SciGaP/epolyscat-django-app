#ifndef PROBLEMLASER_H
#define PROBLEMLASER_H

#include "problemNode.h"
#include "tools.h"

class ProblemLaser : public ProblemNode
{
    std::string _polarization; // z-linear, circular, general
    double _lambda_nm=0; // longest of all wave lenth
    double _intensWcmSq=0;  // peak intensity at longest wave-length
    double _durationOptCyc=0;   // overall duration (by fwhm of individual pulses

public:
    ProblemLaser(ReadInputList& Inp);

    ProblemLaser():ProblemNode("Laser"){}
    std::string definition() const {
        std::vector<std::string> parts=
        {_polarization,tools::str(_lambda_nm,3),tools::str(_intensWcmSq,3),tools::str(_durationOptCyc,3)};
        parts.insert(parts.end(),_inputLines.begin(),_inputLines.end());
        return definitionJoin(parts);
    }
    ProblemLaser(std::string Definition):ProblemLaser(){
        auto parts=definitionSplit(Definition);
        if(parts.size()>3){
            _polarization=parts[0];
            _lambda_nm=tools::string_to_double(parts[1]);
            _intensWcmSq=tools::string_to_double(parts[2]);
            _durationOptCyc=tools::string_to_double(parts[3]);
            _inputLines={parts.begin()+4,parts.end()};
        }
    }

    void addToInput(std::ofstream &Stream) const;

    bool operator==(const ProblemNode & Other) const {
        auto o=dynamic_cast<const ProblemLaser*>(&Other);
        return o and not ((*this)<(*o)) and not ((*o)<(*this));
    };
    bool operator<(const ProblemLaser &Other) const {
        if(_polarization<Other._polarization)return true;
        if(_lambda_nm<Other._lambda_nm)return true;
        if(_intensWcmSq<Other._intensWcmSq)return true;
        if(_durationOptCyc<Other._durationOptCyc)return true;
        return ProblemNode::operator<(Other);
    }
    std::string str() const;
};

#endif // PROBLEMLASER_H
