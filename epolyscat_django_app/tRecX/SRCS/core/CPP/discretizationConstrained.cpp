// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "discretizationConstrained.h"

#include <string>
#include "readInput.h"
#include "operatorMap.h"
#include "algebra.h"
using namespace std;


static void read(ReadInput & Inp, int Line, std::string &AxName, double &LowerEnd, double &UpperEnd, unsigned int &Order){
    Inp.read("Constrain","axis",AxName,"NONE","name of axis to apply constraints",Line);
    Inp.read("Constrain","upperEnd",UpperEnd,Algebra::Infty,"constrain axis to <= upper end",Line);
    Inp.read("Constrain","lowerEnd",LowerEnd,tools::str(-UpperEnd,5),"constrain axis to >= lower end",Line);
    Inp.read("Constrain","order",Order,"1000000","maximal order",Line);
}
bool DiscretizationConstrained::inputs(ReadInput & Inp){
    string axName;
    unsigned int order;
    double lower,upper;
    read(Inp,1,axName,lower,upper,order);
    return axName!="NONE";
}

DiscretizationConstrained::DiscretizationConstrained(const Discretization *D, ReadInput &Inp){

    parent=D;
    name=D->name+"_constrained";
    axis=D->axis;
    continuityLevel=D->continuityLevel;
    string axName;
    unsigned int line=0,order;
    double lower,upper;
    constString="NONE";
    while (line<100){
        line++;
        read(Inp,line,axName,lower,upper,order);
        if(axName=="NONE")break;
        if(order<1)ABORT("in axis "+axName+": cannot constrain to order<1");

        // get axis number
        unsigned int k=0;
        for(;k<axis.size();k++)
            if(axis[k].name==axName)break;
        if(k==axis.size())ABORT("cannot constrain, no such axis: "+axName);
        //        if(axis[k].basDef.size()<2)ABORT("for now, can only constrain finite element axes, is: "+axis[k].name);

        axis[k].constrain(lower,upper,order);

        if(line==1)constString="";
        else constString+=", ";
        constString+=name+"["+tools::str(axis[k].lowerEnd())+","+tools::str(axis[k].upperEnd())+"],order<="+tools::str(order);
    }

    construct();

    _mapFromParent.reset(new OperatorMap(idx(),D->idx()));
    _mapFromParent.reset(new OperatorMap(D->idx(),idx()));
}

