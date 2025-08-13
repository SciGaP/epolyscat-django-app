#ifndef RANDOMPOTENTIAL_H
#define RANDOMPOTENTIAL_H

#include "algebraPiecewise.h"

class ReadInput;

/// random potential on a single axis, expressed as AlgebraPiecewise
class RandomPotential
{
    std::vector<double> _qs;        // boundaries between intervals, including beginning and end (if finite)
    std::vector<std::string> _defs; // definition of potential function on interval
    std::vector<double> _random01;  // numbers in [0,1] that are randomly generated upon construction
    static std::map<std::string,RandomPotential> _list; // list of all random potentials w/o the random numbers resolved

    RandomPotential(std::vector<double>Qs,std::vector<std::string> Defs);
public:
    RandomPotential(){}
    bool operator==(const RandomPotential & Other) const {return _qs==Other._qs and _defs==Other._defs and _random01==Other._random01;}

    /// read random potential from Inp, put into _list by name
    static void read(ReadInput & Inp);

    /// return random potential as Algebra by name (function type algebraFactory)
    static const Algebra* factory(std::string Name);

    static void clearAll(){_list.clear();} ///< removes all potentials from list, names can be re-used
};

#endif // RANDOMPOTENTIAL_H
