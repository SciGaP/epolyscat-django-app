// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "eigenInverseIter.h"

#include "operatorTree.h"
#include "resolvent.h"
#include "coefficients.h"
#include "printOutput.h"

EigenInverseIter::EigenInverseIter(double Epsilon, int Updates, int Power)
    :_epsilon(Epsilon),_updates(Updates),_npower(Power){withSelect("undefined");}

EigenInverseIter &EigenInverseIter::withGuesseigenvalue(std::complex<double> Eguess){
    _select="Nearest[1,"+tools::str(Eguess.real())+","+tools::str(Eguess.imag())+"]";
    return *this;
}


bool EigenInverseIter::converged(const std::vector<std::complex<double> > Eval){
    // either numerically zero change of last two, or less than _epsilon across last three
    if(Eval.size()<2)return false;
    if(std::abs(Eval.back()-Eval[Eval.size()-2])<1.e-14*std::max(1.,std::abs(Eval.back())))return true;
    if(Eval.size()<3)return false;
    for(int k=Eval.size()-3;k<Eval.size();k++)
        if(std::abs(Eval[k-1]-Eval[k])>_epsilon)return false;
    return true;
}

void EigenInverseIter::_compute(){
    if(not _select.find("Nearest[1,"!=0)){
        DEVABORT("need Nearest[1,vReal,vImag], got: "+_select+", call as slv.withGuesseigenvalue(v).compute(Op,Ovr)");
    }

    const OperatorTree* _opTree=dynamic_cast<const OperatorTree*>(_op);
    if(not _opTree)DEVABORT("EigenInverseIter only for OperatorTree");

    if(     not (_op->isComplexSymmetric() and _ovr->isComplexSymmetric()) and
            not (_op->isSelfAdjoint() and _ovr->isSelfAdjoint()))
        PrintOutput::warning("Eigenproblem neither complex symmetric nor selfadjoint - duals will be incorrect");

    Coefficients* in =new Coefficients(_op->iIndex);
    Coefficients* out=new Coefficients(_op->iIndex);
    out->setToRandom();

    std::vector<std::complex<double>> updEvals;
    std::complex<double> ovr=1.,eRes=guessEigenvalue();
    for(int up=0;up<_updates;up++){
        Resolvent res(_opTree,eRes);
        std::vector<std::complex<double>> powEvals(1,eRes);
        for (int p=0;p<_npower;p++){
            _ovr->apply(1.,*out,0,*in);
            res.apply(1./ovr,*in,0.,*out);
            std::complex<double> ovr=_ovr->matrixElement(*out,*out,true);
            std::complex<double> ene=_opTree->matrixElement(*out,*out,true);
            powEvals.push_back(ene/ovr);
            if(converged(powEvals))break;
        }
        updEvals.push_back(powEvals.back());
        if(converged(updEvals))break;
        eRes=updEvals.back();
    }
    if(not converged(updEvals))PrintOutput::warning(Sstr+"inverse iteration not converged, last two values:"
                                                    +updEvals[updEvals.size()-1]+updEvals.back());
    _eigenvalues.clear();
    _eigenvalues.assign(1,updEvals.back());
    _rightVectors.assign(1,out);
    _dualVectors.assign(1,new Coefficients(out->idx()));
    _ovr->apply(1.,*out,0.,*_dualVectors[0]);
    delete in;
}
