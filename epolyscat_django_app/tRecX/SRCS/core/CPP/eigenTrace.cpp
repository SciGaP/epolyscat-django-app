// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "eigenTrace.h"
#include <vector>
#include <complex>
#include <regex>

#include "readInput.h"
#include "asciiFile.h"
#include "parameters.h"
#include "algebra.h"
#include "eigenSolver.h"
#include "operatorTree.h"
#include "operatorDefinition.h"
#include "tools.h"

static std::string funcName="traceF[Q]";

EigenTrace::EigenTrace(ReadInput & Inp)
{
    if(not Inp.found("Trace"))return;

    std::vector<std::string>sVal;
    Inp.texdocuCategoryAdd("Trace","fromQ,toQ,steps,function,eigenvalues",
                           R"tex(
                           Trace eigenvalues of a parameter dependent Hamiltonian.
                           Increment parameter over a range of values, find nearest values by inverse iteration,
                           select value wher wave funstion has most overlap with previous.
                           )tex","90");

    Inp.read("Trace","eigenvalues",sVal,"","guess values, blank-separated list of complex numbers, format 1.2+i3.4 5.6+i7.8 etc.")
            .texdocu(R"tex(
                     Initial guess values for tracing.\\
                     When eigenvectors are degenerate at some \lcode{arg} (see \nameref{docu:Trace:function}), it is useful to determine
                     eigenvalues at slightly different \lcode{arg'} where degeneracy is lifted an start from desired roots there.
                     )tex");
    for(std::string s: sVal)_eStart.push_back(Algebra::constantValue(s));
    std::string fVal;
    Inp.read("Trace","eig-file",fVal,"","guess values on tRecX formated eig-file")
            .texdocu(R"tex(
                     (see \lcode{eigenvalues}))
                     )tex");
    if(fVal!=""){
        if(not folder::exists(fVal))ABORT("Trace:eig-file \""+fVal+"\" does not exists");
        AsciiFile f(fVal);
        std::vector<double> eR,eI;
        f.readCol(eR,0);
        f.readCol(eI,1);
        for(size_t k=0;k<eR.size();k++)_eStart.push_back(std::complex<double>(eR[k],eI[k]));
    };

    std::string funcDef;
    Inp.read("Trace","function",funcDef,"traceF[Q]=(Q)","function for use in operator, func[arg]=(algebra of arg), func[arg] is optional")
            .texdocu(R"tex(
                     Define a function name in the format \lcode{name[arg]=<algebra>} (sec \nameref{sec:class:Algebra}).
                     Strings for name and arg can be chosen freely, choose reasonable memnotics.
                     \lcode{<algebra>} is any admissible algebra string of "arg". If \lcode{name[arg]=} is not specified,
                     "Q" will be interpreted as the argument in \lcode{<algebra>}.
                     )tex");
    std::vector<std::string> part=tools::splitString(funcDef,'=');
    if(part.size()==2){
        funcName=part[0];
        std::string alg=tools::substringReplaceAll(part[1],tools::stringInBetween(part[0],"[","]"),std::string("Q"));
        _parFunc.reset(new Algebra(alg));
    }
    else if(part.size()==1)
        _parFunc.reset(new Algebra(part[0]));
    else
        ABORT("need Trace:function as func[Var]=(algebra of Var) or (algebra of Q), got: "+funcDef);

    if(not _parFunc->isAlgebra())ABORT("not an algebra string: "+funcDef+"\n "+Algebra::failures);
    Inp.read("Trace","fromQ",_parMin,"0.","start value for parameter Q")
            .texdocu(R"tex(
                     First value for "arg" in \nameref{docu:Trace:function}.
                     )tex");
    Inp.read("Trace",  "toQ",_parMax,"0.","last value for parameter Q")
            .texdocu(R"tex(
                     Last value for "arg" in \nameref{docu:Trace:function}.
                     )tex")
            .inputUnits(_parUnits);

    Inp.read("Trace","steps",_nSteps,_parMin<_parMax?"11":"1","number of steps in [from,to]")
            .texdocu(R"tex(
                     Number of equidistant steps for "arg" in the interval [\nameref{docu:Trace:fromQ},\nameref{docu:Trace:toQ}].
                     )tex");
    if(_eStart.size()){
        if(_parMin==_parMax)PrintOutput::warning("zero parameter range, will only calculate for first guess");
        Parameters::addResetable(funcName,_parFunc->val(_parMin));
        Algebra::addUpdatableConstant(funcName,_parFunc->val(_parMin)); // let be updated, if needed
    }
    Inp.read("Trace","ovrErr",_ovrErr,"Infty","stop when overlaps error of subsequent eigenvectors exceeds ovrErr")
            .texdocu(R"tex(
                     Stop tracing, if overlap of subsequent wave functions $|\l \Psi(\text{arg}_i)|\Psi(\text{arg}_{i+1}\r|>$\lcode{overErr}
                     This can (and will) happen near (avoided-) crossings of roots.
                     )tex");
    print();
}


void EigenTrace::print() const{
    std::string var=tools::stringInBetween(funcName,"[","]");
    PrintOutput::title("Trace eigenvalues");
    PrintOutput::paragraph();
    PrintOutput::newRow();
    PrintOutput::rowItem("name");
    PrintOutput::rowItem("function(Q="+var+")");
    PrintOutput::rowItem("from "+var);
    PrintOutput::rowItem("to "+var);
    PrintOutput::newRow();
    PrintOutput::rowItem(funcName);
    PrintOutput::rowItem(_parFunc->definition());
    PrintOutput::rowItem(_parMin);
    PrintOutput::rowItem(_parMax);
    PrintOutput::paragraph();
    if(_eStart.size()<6)
        PrintOutput::lineItem("initial guess(es)",tools::str(_eStart));
    else{
        double eMin=DBL_MAX,eMax=-DBL_MAX;
        for(auto e: _eStart)eMin=std::min(eMin,e.real());
        for(auto e: _eStart)eMax=std::max(eMax,e.real());
        PrintOutput::lineItem("initial guesses",Str("","")+_eStart.size()+" values in ["+eMin+","+eMax+"]");
    }
}

bool EigenTrace::run(std::string OpDef, const Index *Idx, std::string RunDir){
    if(not _eStart.size())return false;

    // create operator
    if(OpDef.find(funcName)==std::string::npos)
        PrintOutput::warning("parameter \""+funcName+"\" not found in\n"+OpDef);
    _op.reset(new OperatorTree("trace",OperatorDefinition(OpDef,Idx->hierarchy()).str(),Idx,Idx));

    AsciiFile fR(RunDir+"traceEigReal"),fI(RunDir+"traceEigImag");
    std::vector<std::vector<std::complex<double>>> traces(1);
    for(auto e: _eStart){
        std::string filnam=RunDir+"traceEig";
        if(_eStart.size()>1)filnam+=tools::str(e.real(),4);
        AsciiFile traceFile(filnam);
        std::vector<std::string>comm;
        comm.push_back("trace eigenvectors for parameter function "+_parFunc->definition());
        std::string par=tools::stringInBetween(funcName,"[","]");
        if(_parUnits!="DEFAULT_SYSTEM")par+="("+_parUnits+")";
        comm.push_back("   Parameter     "+par+"         Re(E)            Im(E)        Overlap w. Previous");
        traceFile.writeComments(comm);
        traces.push_back(std::vector<std::complex<double>>());
        trace(e,traces[0],traces.back(),traceFile);

        std::vector<std::vector<double>>cols;
        std::vector<std::string>comms;
        traceFile.readComments(comms);
        traceFile.readCols(cols,{0,1,2,3},",");
        if(filnam!="traceEig"){
            if(folder::exists(RunDir+"traceEigReal")){
//                AsciiFile fR(Dir+"traceEigReal"),fI(Dir+"traceEigImag");
                if(not fR.addCols({cols[2]},{}))
                    PrintOutput::warning("first column on traceEig does not match first on "+filnam+" -- not added");
                for(auto &a: cols[3])a=-a;
                fI.addCols({cols[3]},{});
            }
            else {
                comms.back()="# "+par+"     Re(E)";
                fR.writeComments(comms);
                fR.writeCols({cols[0],cols[2]});
                comms.back()="# "+par+"    -Im(E)";
                fI.writeComments(comms);
                for(auto &a: cols[3])a=-a;
                fI.writeCols({cols[0],cols[3]});
            }
        }
    }
    PrintOutput::title("all traces on (real,imag)= "+fR.name()+", "+fI.name());
    return true;
}

void EigenTrace::trace(std::complex<double> Eguess, std::vector<std::complex<double>> & Par, std::vector<std::complex<double>> & Eval,AsciiFile & File){
    Coefficients prev(_op->iIndex,0.);
    PrintOutput::title(Sstr+"Trace root starting from "+Eguess);
    PrintOutput::paragraph();
    PrintOutput::newRow();
    std::string par=tools::stringInBetween(funcName,"[","]");
    if(_parUnits!="DEFAULT_SYSTEM")par+="("+_parUnits+")";
    PrintOutput::rowItem("   "+par);
    PrintOutput::rowItem("    "+funcName+"  ");
    PrintOutput::rowItem("      Re(E)   ");
    PrintOutput::rowItem("      Im(E)   ");
    PrintOutput::rowItem("Overlap w. previous");
    for(int k=0;k<_nSteps;k++){
        EigenSolver slv;
        // first step is initial - get only single value
        slv.withSelect("Nearest["+tools::str(k?1:1)+","+tools::str(Eguess.real())+","+tools::str(Eguess.imag())+"]");

        //  updata trace parameter
        double traPar=_parMin+k*(_parMax-_parMin)/double(std::max(_nSteps-1,1));
        double couple=_parFunc->val(traPar).real();
        Par.push_back(std::complex<double>(traPar,couple));
        Algebra::addUpdatableConstant(funcName,couple); // update elementary algebra
        Parameters::update();

        PrintOutput::outputLevel("low");
        slv.compute(_op.get());
        slv.normalize();
        PrintOutput::outputLevel("restore");

        // select best match
        std::vector<std::complex<double>> ovr;
        for(size_t l=0;l<slv.eigenvalues().size();l++)
            ovr.push_back(_op->iIndex->overlap()->matrixElement(prev,*slv.rightVectors()[l]));


        int kMatch=0;
        for(size_t l=0;l<ovr.size();l++)if(std::abs(ovr[l])>std::abs(ovr[kMatch]))kMatch=l;
        if(k>0 and std::abs(std::abs(ovr[kMatch])-1.)>_ovrErr){
            PrintOutput::warning(Sstr+"lost overlap to previous step - use denser steps, last overlap ="+ovr);
            break;
        }

        prev=*slv.rightVectors()[kMatch];
        Eval.push_back(slv.eigenvalues()[kMatch]);
        if(Eval.size()>1)Eguess=2.*Eval.back()-Eval[Eval.size()-2];
        else Eguess=Eval.back(); // guess value for next step
        Par.back()=std::complex<double>(Units::convert(Par.back().real(),"DEFAULT_SYSTEM",_parUnits),Par.back().imag());
        PrintOutput::newRow();
        PrintOutput::rowItem(Par.back().real(),5);
        PrintOutput::rowItem(Par.back().imag(),8);
        PrintOutput::rowItem(Eval.back().real(),10);
        PrintOutput::rowItem(Eval.back().imag(),4);
        PrintOutput::rowItem(std::abs(ovr[kMatch]),4);
        PrintOutput::flush();
        std::vector<double> row={Par.back().real(),Par.back().imag(),Eval.back().real(),Eval.back().imag(),std::abs(ovr[kMatch])};
        File.writeRow({Par.back().real(),Par.back().imag(),Eval.back().real(),Eval.back().imag(),std::abs(ovr[kMatch])});
    }
    PrintOutput::paragraph();
    PrintOutput::message("eigenvalue trace on "+File.name());
}
