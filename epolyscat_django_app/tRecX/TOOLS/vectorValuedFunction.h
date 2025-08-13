// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef VECTORVALUEDFUNCTION_H
#define VECTORVALUEDFUNCTION_H

#include <vector>
#include <complex>
#include <map>
#include <memory>
#include "abort.h"

class VectorValuedFunction
{
    std::string _funcName;
    static std::map<std::string,std::shared_ptr<const VectorValuedFunction> > _list;
    std::vector<double> _lowLim;
    std::vector<double> _upLim;
public:
    VectorValuedFunction(){}
    virtual std::vector<std::complex<double> > operator()(std::vector<double> X) const=0;
    virtual std::string coordinates() const=0; ///< coordinates, e.g. X.Y.Z or Phi.Eta.Rn
    virtual unsigned int length() const=0; ///< length of vector

    virtual std::string info(int Number=-1) const;

    /// add to list by name
    static void add(std::string Name,std::shared_ptr<const VectorValuedFunction> Func){
        if(_list.count(Name) and _list[Name].get()!=Func.get())ABORT("multiple definition of "+Name);
        _list[Name]=Func;
    }

    /// fetch from list by name
    static const VectorValuedFunction* get(std::string Name){
        auto p=_list.find(Name);
        if(p!=_list.end())return p->second.get();
        return 0;
    }

    std::string name() const{return _funcName;}
    virtual bool validArguments(std::vector<double> X, std::string Message="") const;
    void abortInvalidArguments(std::string Message="") const {
        ABORT("arguments invalid for"+name()+" (coordinates "+coordinates()+"): "+Message);
    }
    static std::string list(); ///< lists of all currently available function names
    static void print(); /// brief info about available functions
};

#endif // VECTORVALUEDFUNCTION_H
