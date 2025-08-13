// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef LABELLED_H
#define LABELLED_H

#include <string>
#include <limits.h>
#include "tools.h"
#include "mpiWrapper.h"

/// labels for class T: label will be created only upon first call to label()
template <typename T>
class Labelled
{
    static long long _nextLabel;
    static std::map<const Labelled<T>*,long> _listOfLabels; // memory-light version
public:
    ~Labelled(){_listOfLabels.erase(this);}

    /// compose a hash string from numerical label and possible overflow
    std::string hash() const {return "#" + tools::str(label(false));}


    /// return label, check for possible overflow
    int label(bool All=false) const {
        // altough long SHOULD be default constructed to =0, we do not rely on it
        auto p=_listOfLabels.find(this);
        if(p==_listOfLabels.end()){
            if(++_nextLabel>=LONG_MAX)DEVABORT("Label overflow: Change implementation.");
            _listOfLabels[this]=_nextLabel;
            p=_listOfLabels.find(this);
        }
        if(p->second > INT_MAX)DEVABORT("Label overflow: Use hash.");
        return int(p->second);
    }
};

template<typename T>
long long int Labelled<T>::_nextLabel=1;

template<typename T>
std::map<const Labelled<T>*,long> Labelled<T>::_listOfLabels;

#endif // LABELLED_H
