// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisMatMatrix.h"

#include "readInput.h"
#include "printOutput.h"
#include "basisMatOperator.h"
#include "basisMatNumbers.h"
#include "algebra.h"
#include "basisAbstract.h"
#include "basisMat1D.h"
#include "basisSub.h"
#ifdef _USE_HACC_
#include "basisMatCI.h"
#endif


std::map<std::string,std::shared_ptr<BasisMatMatrix> > BasisMatMatrix::_list;

void BasisMatMatrix::add(std::string Def,const Eigen::MatrixXcd & Mat){
    if(_list.count(Def)){
        if(not (_list[Def]->_mat==Mat))
            ABORT("matrix "+Def+" defined previously with different values:\n"
                  +EigenTools::str(_list[Def]->_mat)
                  +"\npresent entries\n"
                  +EigenTools::str(Mat));
    }
    else{
        _list[Def].reset(new BasisMatMatrix());
        _list[Def]->_mat=Mat;
    }
}

void BasisMatMatrix::read(ReadInput & Inp){
    std::string namePrev="NONE",name,kind;

    // advance to Number'th definition
    int line=0;
    while(not Inp.endCategory("Matrix",++line)){
        Inp.read("Matrix","name",name,"","free-chosen name of matrix",line);
        if(name!=namePrev){
            Inp.read("Matrix","kind",kind,"matrixElements","Input mode:"
                     +std::string(" matrixElements...elemente-wise input")
                     +std::string(", diagonal...blank-separated list of diagonal values")
                     +std::string(", operator...operator for given basis")
                     ,line);

            if(kind.find("matrixElements")==0)_list[name].reset(new BasisMatNumbers(Inp,line));
            else if(kind.find("diagonal")==0)_list[name].reset(new BasisMatNumbers(Inp,line));
            else if(kind.find("operator")==0)_list[name].reset(new BasisMatOperator(kind));
            else ABORT("undefined Matrix: kind="+kind);
        }
    }
}


const BasisMatMatrix* BasisMatMatrix::factory(std::string Op, const BasisAbstract* IBas, const BasisAbstract * JBas)
{
    if(Op.find("<")==std::string::npos)ABORT("operator definition must contain <...>, is: "+Op);
    std::complex<double> preFac=preFactor(Op);

    std::string Op0=tools::stringInBetween(Op,"<",">");
    std::string OpWithBasis=Op+":"+IBas->str()+"|"+JBas->str();

    BasisMatMatrix * m=0;

    const BasisAbstract* iSup=BasisSub::superBas(IBas);
    const BasisAbstract* jSup=BasisSub::superBas(JBas);

    if(iSup->hybrid() or jSup->hybrid()){
        if(IBas!=JBas)DEVABORT("need equal left and right hybrid basis, got:"+IBas->str()+" "+JBas->str());
        BasisMat1D m1d(Op,iSup->size(),jSup->size());
        Eigen::MatrixXcd mm;
        if(m1d.isEmpty()){
            OpWithBasis="<allOnes>:"+tools::str(IBas->size())+"|"+tools::str(JBas->size());
            mm=Eigen::MatrixXcd::Constant(IBas->size(),JBas->size(),1.);
        }
        else
            mm=m1d.mat();
        if(not _list.count(OpWithBasis)){
            _list[OpWithBasis].reset(new BasisMatMatrix());
            _list[OpWithBasis]->_mat=mm;
        }
        m=_list[OpWithBasis].get();
    }

    else if(Op=="<allOnes>"){
        std::string onesSize=Op+"<allOnes>:"+tools::str(IBas->size())+"|"+tools::str(JBas->size());
        if(not _list.count(onesSize)){
            _list[onesSize].reset(new BasisMatMatrix());
            _list[onesSize]->_mat=Eigen::MatrixXcd::Constant(IBas->size(),JBas->size(),1.);
        }
        m=_list[onesSize].get();
    }

#ifdef _USE_HACC_
    else if (iSup->ci() or iSup->ci()){
            BasisMatCI mci(Op,IBas,JBas);
            _list[OpWithBasis].reset((m=new BasisMatNumbers(mci.mat())));
    }
#endif

    else if(_list.count(OpWithBasis)){
        m=_list[OpWithBasis].get();
    }else if(_list.count(Op)){
        m=_list[Op].get();
    }else if(_list.count(Op0)){
        m=_list[Op0].get();
    } else {
        BasisMat1D m1d(Op,IBas,JBas);
        // save for re-use, but not for random-type operators
        if(not m1d.isEmpty() and Op.find("random")==std::string::npos)
            _list[OpWithBasis].reset((m=new BasisMatNumbers(m1d._mat)));
    }
    if(not m)return 0;

    if(m->mat().size()==0 and dynamic_cast<BasisMatOperator*>(m)){
        dynamic_cast<BasisMatOperator*>(m)->setup(IBas,JBas);
        _list[OpWithBasis].reset(new BasisMatNumbers(m->_mat));
    }

    if(preFac!=1.)_list[OpWithBasis].reset(new BasisMatNumbers(m->_mat*preFac));

    if(m->_mat.size()==0)return 0;

    // check dimensions
    if(m->_mat.rows()<IBas->size() or m->_mat.cols()<JBas->size())
        ABORT(Str("matrix dimensions for"," ")+Op+"too small:"
              +m->_mat.rows()+"x"+m->_mat.cols()+"<"+IBas->size()+"x"+JBas->size()+IBas->strDefinition()+JBas->str());
    if(m->_mat.rows()!=IBas->size() or m->_mat.cols()!=JBas->size()){
        PrintOutput::warning(
                    Str("matrix dimensions for"," ")+Op+"exceed basis sizes (using lowest sub-block):"
                    +m->_mat.rows()+"x"+m->_mat.cols()+">"+IBas->size()+"x"+JBas->size());
        m->_mat.conservativeResize(IBas->size(),JBas->size());
    }
    return m;
}


const UseMatrix BasisMatMatrix::useMat() const {
    return UseMatrix::UseMap(const_cast<BasisMatMatrix*>(this)->_mat.data(),_mat.rows(),_mat.cols());
}

std::complex<double> BasisMatMatrix::preFactor(std::string Op){
    std::string fac=Op.substr(0,Op.find("<"));
    if(fac=="")return 1.;
    if(fac=="+")return 1.;
    if(fac=="-")return -1.;
    Algebra facAlg(fac);
    if(not facAlg.isAlgebra())ABORT("prefactor in "+Op+" is not a valid Algebra");
    if(not facAlg.isAlgebraOfConsts())ABORT("only Algebra of constants as prefactor, found: "+fac+" in "+Op);
    return facAlg.val(0.);
}
