// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "tools.h"
#include "tools.h"
#include <deque>
#include "functionOneArg.h"
#include "algebra.h"

#ifdef _WIN32
// for alternative tokens for windows ( e.g. and, or, not, ... ) 
#include <iso646.h>
#endif

class ReadInput;

/** \ingroup OperatorData */

//! \brief interface for Parameter's or functions that can be updated
//!
//! Parameters has a table of (pointers to) Updatable's. Whenever Parameters::update is invoked,
//! Updatable::update will be invoked on every element.
class Updatable
{
public:
    virtual ~Updatable()=0;
    virtual void update(double time)=0;
};


/** \ingroup OperatorData */

/// parameters that are use in operator definitions
class Parameters {
//    friend class Operator;
    friend class OperatorDefinition;
public:

    /// parameter function: used to update parameter for given time
    typedef std::complex<double> (*parFunction)(const double time);

protected:
    std::string name;                     ///< unique name for parameter
    std::complex<double> plusValue;       ///< positive current parameter value
    std::complex<double> minusValue;      ///< negative current parameter value
    parFunction updateFunction;           ///< OBSOLESCENTE how to compute given parameter for given time
    const FunctionOneArg * updateFunctionArg;   ///< OBSOLESCENT alternate, new form of function update
    const Algebra* algebra;              ///< compute parameter as algebra
    bool resetable;                        ///< indicates that parameter can be changed

    // NOTE: must use deque, as vector moves its data and pointers loose meaning
    static std::deque<Parameters> table; //!< list of all factors that can be pointed to by floors
    static double lastUpdateTime;           //!< time used in the last update
    static std::deque<Updatable*> updatables; //!< table of pointers to all objects that will (and therefore can) be updated
    static Parameters noParameter;
    static bool parsForcedTo1;

public:
    Parameters():name("NOT_SET"),updateFunction(0),updateFunctionArg(0),algebra(0){}

    //! arguments: unique string (e.g."Apot", "z-Field"), current value, function pointer for updating current value
    Parameters(std::string Name, std::complex<double> CurrentValue,parFunction Update,bool Reset)
        :name(Name),plusValue(CurrentValue),minusValue(-CurrentValue),updateFunction(Update),updateFunctionArg(0),algebra(0),resetable(Reset){}

    //! arguments: unique string to identify FunctionOneArg
    Parameters(std::string Name):name(Name),plusValue(1.),minusValue(-1.),updateFunction(0),
        updateFunctionArg(FunctionOneArg::get(Name)),resetable(false){}

    /// set default pre-defined factors
    static void defaults();

    /// list of parameters that are functions of parameters
    static void setSpecial();
    static void updateSpecial();

    /// add to list of parameters: string, current value, update function (or null-pointer for constant parameters)
    /// NOTE: from classes, only static class functions can be added (NOT member functions)
    static void add(const Algebra* Alg);
    static void add(std::string Name, const std::complex<double> CurrentValue, parFunction Update=0);
    static void add(std::string Name);
    static void addResetable(std::string Name,const std::complex<double> CurrentValue=std::complex<double>(1.,0.));
    static void reset(std::string Name,const std::complex<double> CurrentValue);
    static void read(std::string Cat,std::string Name,std::string Default, ReadInput & in); ///< read constant parameter directly from input

    static void add(Updatable* ptUpdatable); //!< add an Updatable to the table

    static double currentTime() {return lastUpdateTime;}
    static void update(const double time, bool Special=true); ///< (re-)evaluate all time-dependent parameters for given time
    static void update();  ///< update if isAlgebraOfUpdatables()
    static void updateToOne();             ///< set all time-dependent parameters to 1 (use with discretion, restore() when no longer needed)
    static void restoreToTime(){parsForcedTo1=false;update(lastUpdateTime);} ///< update parameters to the previously set time (cleanup after updateToOne())

    /// return parameter from table
    static Parameters & tableParameter(const std::string Name, bool Abort=true);

    /// true if value can be updated
    static bool isFunction(std::string Name);
    static bool isDefined(std::string Name);

    /// return pointer to parameter value
    static std::complex<double> *pointer(std::string Name);
    static void show();///< list available parameters
};

#endif
