// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "eigenTools.h"
#include "eigenSolverAbstract.h"
#include "coefficients.h"
#include "operatorTree.h"
#include "operatorSubspace.h"

#include "tools.h"
#include "printOutput.h"
#include "tRecXchecks.h"
#include "algebra.h"
#include "readInput.h"
#include "operatorDefinition.h"
#include "operatorSubspace.h"
#include "projectSubspace.h"

//DEBUG
#include "mpiWrapper.h"
#include "parallelOperator.h"

using namespace std;
std::string EigenSolverAbstract::method;
std::string EigenSolverAbstract::defaultMethod;

void EigenSolverAbstract::makeSelect(int Nev, double Emin, double Emax, bool ExcludeRange) {

    std::string sMin=tools::str(Emin,5,DBL_MAX/2);
    std::string sMax=tools::str(Emax,5,DBL_MAX/2);
    std::string sNev=Nev>INT_MAX/2?Algebra::Infty:tools::str(Nev);

    if(sMin.find("-"+Algebra::Infty)==0 and  sMin.find(Algebra::Infty)==0 and sNev==Algebra::Infty)
        _select="All";
    else if(sNev==Algebra::Infty){
        _select="Range["+sMin+","+sMax+"]";
        if(ExcludeRange)_select="Outside"+_select;
    }
    else if(sMin!="-"+Algebra::Infty or sMax!=Algebra::Infty){
        _select="Range["+sMin+","+sMax+","+sNev+"]";
        if(ExcludeRange)_select="Outside"+_select;
    }
    else if(sMin=="-"+Algebra::Infty)
        _select=ExcludeRange?"LargeAbs["+sNev+"]":"SmallReal["+sNev+"]";
    else if(sMax==Algebra::Infty)
        _select=ExcludeRange?"SmallReal["+sNev+"]":"LargeAbs["+sNev+"]";
    else
        ABORT("undefined combination of parameters for spectral selection: Nev,Emin,Emax = "+sNev+","+sMin+","+sMax);
}


std::complex<double> EigenSolverAbstract::guessEigenvalue() const{
    if(_select.find("Nearest[")==0){
        std::vector<std::string> part=tools::splitString(tools::stringInBetween(_select,"[","]"),',');
        return std::complex<double>(tools::string_to_double(part[1]),tools::string_to_double(part[2]));
    }
    else DEVABORT("no guessEigenvalue for "+_select);
    return 0.;
}

void EigenSolverAbstract::readControls(ReadInput &Inp){
    Inp.obsolete("Eigensolver","method","use Eigen: method instead");
    Inp.texdocuCategoryAdd("Eigen","select,method,vectors,shif","control which eigenvalues are computed, and how","00,02,05,22,90");
    Inp.read("Eigen","method",defaultMethod,"auto","options: auto...(default) reasonable choice,Lapack,Arpack,NonLin",1,"eigenMethod")
            .texdocu(R"tex(
                     Not all methods will work with all \nameref{docu:Eigen:select}'s.
                     \begin{itemize}
                     \item[auto] selects a method depending on Hamiltonian structure and size
                     \item[Lapack] a full matrix will be formed and standard Lapack will be used for solving.
                     recommended up to problem sizes $\lesssim 3000$ for the diagonal blocks of the Hamiltonian matrix.
                     \item[Arpack] Lanczos Arnoldi from the Arpack package.
                     \item[NonLin] [EXPERT] for non-linear eigenproblems
                     \end{itemize}
                     )tex");
    method=defaultMethod;

}
void EigenSolverAbstract::setMethod(std::string Method){
    if(method!=defaultMethod)PrintOutput::DEVwarning("setting new EigenSolver method="+Method+" but previous != default");
    method=Method;
}

std::string EigenSolverAbstract::autoMethod(const string Method){
    if(Method!="auto")return Method;
    if(_op==0)return Method;
    int size=_op->iIndex->sizeCompute();
    if(size<2000 or numberEigenv()*5>size)return "Lapack";
    if(_select.find("Nearest[1,")==0)return "InverseIteration";
    if(_select.find("Nearest[")==0)return "ArpackShiftInverse";
    return "Arpack";
}

EigenSolverAbstract::~EigenSolverAbstract(){
    Coefficients::DEBUGstatus="SLV";
    for (auto v: _rightVectors)delete v;
    for (auto v: _dualVectors)delete v;
    for (auto v: _leftVectors)delete v;
    Coefficients::DEBUGstatus="";

}

EigenSolverAbstract::EigenSolverAbstract(const OperatorAbstract *Op):
    _serial(true),
    _op(Op),
    _ovr(0),
    _computeLeftVectors(false),
    _computeRightVectors(true),
    _sort("SmallReal"),
    _epsVerified(DBL_MAX)
{
    makeSelect(Op?int(Op->iIndex->sizeStored()):INT_MAX,-DBL_MAX,DBL_MAX,false);
}

EigenSolverAbstract & EigenSolverAbstract::compute(const OperatorAbstract *Op, const OperatorAbstract *Ovr){
    //if(_op!=0)ABORT("Operator alread set - call copute w/o arguments");
    LOG_PUSH("compA1");
    _op=Op;

    if(_ovr==0){
        if(Ovr==0)_ovr=Op->iIndex->overlap();
        else _ovr=Ovr;
    }
    else if(Ovr!=0)
        //ABORT("Overlap alread set - call compute w/o Ovr argument");
        if(_op->iIndex!=_ovr->iIndex)ABORT("bad");

    _rightVectors.clear();
    _leftVectors.clear();
    _dualVectors.clear();
    _eigenvalues.clear();

    _epsVerified=DBL_MAX;
    _compute();
    LOG_POP();
    LOG_PUSH("compA2");
    _postProcess();
    LOG_POP();
    LOG_PUSH("compA3");
    normalize();
    LOG_POP();
    return *this;
}

void EigenSolverAbstract::_postProcess(){
    if(dynamic_cast<const OperatorSubspace*>(_op)){
        // remove possible sub-space content from rightVecors
        rightVectors(); // ensure existence of vectors on this level
        dualVectors();
        const OperatorSubspace* os=dynamic_cast<const OperatorSubspace*>(_op);
        for (auto &v: _rightVectors)os->projector()->apply(-1,*v,1,*v);
        for (auto &v: _dualVectors)os->projector()->applyDual(-1,*v,1,*v);
    }
}

std::map<std::string,std::string> selectPars(std::string Select){
    std::vector<std::string> p=tools::splitString(tools::stringInBetween(Select,"[","]"),',');
    std::map<std::string,std::string> m;
    m["eMin"]="-"+Algebra::Infty;
    m["eMax"]=Algebra::Infty;
    m["nEv"]=Algebra::Infty;
    m["Exclude"]="0";
    if(0==Select.find("All")){
    }
    else if(0==Select.find("Range[")){
        if(p.size()<2)ABORT("need Range[Emin,Emax{,Nmax}], got: "+Select);
        m["eMin"]=p[0];
        m["eMax"]=p[1];
        if(p.size()>2)m["nEv"]=p[2];
    }
    else if(0==Select.find("SmallReal[")){
        if(p.size()!=1)ABORT("need SmallReal[N], got: "+Select);
        m["nEv"]=p[0];
    }
    else if(0==Select.find("LargeAbs[")){
        if(p.size()<1)ABORT("need LargeAbs[N], got: "+Select);
        m["nEv"]=p[0];
    }
    else if(0==Select.find("OutsideRange[")){
        if(p.size()!=2)ABORT("need OutsideRange[Emin,Emax], got: "+Select);
        m["eMin"]=p[0];
        m["eMax"]=p[1];
        m["Exclude"]="1";
    }
    else if(0==Select.find("Nearest[")){
        if(p.size()!=3)ABORT("need Nearest[N,Ereal,Eimag], got: "+Select);
        m["nEv"]=p[0];
        m["Ereal"]=p[1];
        m["Eimag"]=p[2];
    }
    else if(0==Select.find("Rectangle[")){
        if(p.size()!=4)ABORT("need Rectangle[realMin,realMax,imagMin,imagMax], got: "+Select);
        m["eMin"]=p[0];
        m["eMax"]=p[1];
        m["iMin"]=p[2];
        m["iMax"]=p[3];
    }
    else
        ABORT("cannot interprete eigenvalue selection "+Select);
    return m;
}

double EigenSolverAbstract::eMin() const {return Algebra(selectPars(_select)["eMin"]).val(0.).real();}
double EigenSolverAbstract::eMax() const {return Algebra(selectPars(_select)["eMax"]).val(0.).real();}
int EigenSolverAbstract::numberEigenv() const {std::string sN=selectPars(_select)["nEv"];
                                               return sN==Algebra::Infty?INT_MAX:int(Algebra(selectPars(_select)["nEv"]).val(0.).real());}
bool EigenSolverAbstract::excludeRange() const {return Algebra(selectPars(_select)["nEv"]).val(0.).real()==0;}

bool EigenSolverAbstract::isSelfAdoint() const {
    return _op->isSelfAdjoint() and _op->idx()->overlap()->isSelfAdjoint();
}
bool EigenSolverAbstract::isComplexSymmetric() const {
    return _op->isComplexSymmetric() and _op->idx()->overlap()->isComplexSymmetric();
}


void EigenSolverAbstract::withSelect(std::string Select){_select=Select;}
std::vector<std::complex<double> >EigenSolverAbstract:: eigenvalues() { return _eigenvalues; }
std::vector<Coefficients* > EigenSolverAbstract::leftVectors() { return _leftVectors; } ///< left hand eigenvectors: (left_i)^* A = e_i (left_i)^* S
std::vector<Coefficients* > EigenSolverAbstract::rightVectors() { return _rightVectors; } ///< right hand eigenvectors Op * right_j = S * right_j eval_j
std::vector<Coefficients* > EigenSolverAbstract::dualVectors() {return _dualVectors;} ///< Dual vectors: Dual[m].innerProduct(Evec[n])=delta[m,n]

std::vector<Coefficients* > EigenSolverAbstract::symmetricDuals() {
    if(_dualVectors.size()!=_rightVectors.size()){
        if(_dualVectors.size()!=0)DEVABORT("number of dual vectors does not match eigenvectors");
        for(const Coefficients* c: _rightVectors){
            _dualVectors.push_back(new Coefficients(_op->idx()));
            _op->idx()->overlap()->apply(1.,*c,0.,*_dualVectors.back());
        }
        _dualVectors.back()->makeContinuous();
        if(isSelfAdoint())
            _dualVectors.back()->conjugate();
        else if (not isComplexSymmetric())
            PrintOutput::warning("problem is neither self adjoint nor complex symmetric - duals and projector incorrect");

    }
    orthonormalizeDegenerate(1.e-12,_eigenvalues,_dualVectors,_rightVectors);

    return _dualVectors;
}

std::string EigenSolverAbstract::orthonormalize(std::vector<Coefficients*>&Dual,std::vector<Coefficients*>&Evec,double EliminateSingular){

    if(Dual.size()!=Evec.size())ABORT("numbers dual eigenvectors not equal right eigenvectors");
    std::string mess("");
    for(size_t i=0;i<Dual.size();i++){
        for(size_t j=0;j<i;j++){
            complex<double> norm = sqrt(Dual[j]->innerProduct(Evec[j], true));
            if(abs(norm)<1.e-12 and EliminateSingular<0.)
                ABORT(Sstr+"cannot orthonormalize - vectors are (near-)linear dependent, pair:"+i+j);
            if(abs(norm)<EliminateSingular){
                norm=1.;
                Dual[j]->setToZero();
                Evec[j]->setToZero();
            }
            complex<double>oij=Dual[j]->innerProduct(Evec[i],true) / norm;
            Evec[i]->axpy(-oij,Evec[j]);
            oij=Dual[i]->innerProduct(Evec[j],true) / norm;
            Dual[i]->axpy(-oij,Dual[j]);
        }
        complex<double>norm=sqrt(Dual[i]->innerProduct(Evec[i],true));
        if(abs(norm)<1.e-12 and EliminateSingular<0.)
            ABORT(Sstr+"cannot orthonormalize - vectors are (near-)linear dependent, zero at:"+i);
        if(abs(norm)<EliminateSingular){
            norm=1.;
            Dual[i]->setToZero();
            Evec[i]->setToZero();
        }
        else {
            std::complex<double> phaseAtMax=Evec[i]->cMaxNorm();
            phaseAtMax/=std::abs(phaseAtMax);
            Evec[i]->scale(std::conj(phaseAtMax)/std::abs(norm));
            norm=Dual[i]->innerProduct(Evec[i],true);
            Dual[i]->scale(1./norm);
        }
    }
    for(size_t i=Evec.size();i>0;i--){
        if(Evec[i-1]->isZero()){
            delete Evec[i-1];
            delete Dual[i-1];
            Evec.erase(Evec.begin()+1-i);
            Dual.erase(Dual.begin()+1-i);
            mess=" "+tools::str(i)+mess;
        }
    }
    return mess==""?"":mess.substr(1);
}


void EigenSolverAbstract::orthonormalizeDegenerate(double Eps,const vector<complex<double> > & Eval,
                                                   const vector<Coefficients*> & Dual,
                                                   const vector<Coefficients*> & Rvec)
{
    // orthogonalize on near-degenerate subspaces, generate duals in the process
    vector<bool> done(Eval.size(),false);
    for(size_t k=0;k<Eval.size();k++){
        if(done[k])continue;
        vector<Coefficients*> dual,rvec;
        for(size_t l=k;l<Eval.size();l++){
            if(done[l])continue;
            if(abs(Eval[k]-Eval[l])<Eps){
                done[l]=true;
                rvec.push_back(Rvec[l]);
                dual.push_back(Dual[l]);
            }
        }
        orthonormalize(dual,rvec);
    }

}

void EigenSolverAbstract::select(string Kind){
    LOG_PUSH("select");
    LOG_PUSH("select1");
    eigenvalues();
    rightVectors();
    leftVectors();
    dualVectors();
    LOG_POP();

    if(_eigenvalues.size()==0)
        PrintOutput::DEVwarning("no eigenvalues to select from - may need to call eigenvalues() before select");
    int eSize=_eigenvalues.size();
    std::string mess;
    if(Kind=="All"){
        sort("SmallReal");
    }

    else if(Kind.find("SmallReal[")==0){
        LOG_PUSH("select2");
        size_t nend=Algebra(tools::stringInBetween(Kind,"[","]")).val(0.).real();
        sort("SmallReal");
        if(nend<_eigenvalues.size()){
            for(size_t  k=nend;k<_rightVectors.size();k++)delete _rightVectors[k];
            for(size_t  k=nend;k<_leftVectors.size();k++)delete _leftVectors[k];
            for(size_t  k=nend;k<_dualVectors.size();k++)delete _dualVectors[k];
            _rightVectors.resize(nend);
            _leftVectors.resize(nend);
            _dualVectors.resize(nend);
            _eigenvalues.resize(nend);
            mess=Sstr+"- max energy"+_eigenvalues.back();
        }
        LOG_POP();
    }

    else if(Kind.find("Nearest[")==0){
        vector<string> rect=tools::splitString(tools::stringInBetween(Kind,"[","]"),',');
        if(rect.size()!=3)ABORT("need Nearest[N,Ereal,Eimag], got: "+Kind);
        size_t nend=Algebra(rect[0]).val(0.).real();
        std::complex<double> eguess(Algebra(rect[1]).val(0.).real(),Algebra(rect[2]).val(0.).real());
        for(auto &e: _eigenvalues)e-=eguess;
        sort("SmallAbs");
        for(auto &e: _eigenvalues)e+=eguess;
        if(nend>=_eigenvalues.size())return;
        for(size_t k=nend;k<_rightVectors.size();k++)delete _rightVectors[k];
        for(size_t  k=nend;k<_leftVectors.size();k++)delete _leftVectors[k];
        for(size_t  k=nend;k<_dualVectors.size();k++)delete _dualVectors[k];
        _rightVectors.resize(nend);
        _leftVectors.resize(nend);
        _dualVectors.resize(nend);
        _eigenvalues.resize(nend);
        mess=Sstr+"- largest distance for"+_eigenvalues.back();
        sort("SmallReal");
    }

    else if(Kind.find("Rectangle[")==0 or Kind.find("Range[")==0){
        vector<string> rect=tools::splitString(tools::stringInBetween(Kind,"[","]"),',');
        vector<double>lim;
        if(Kind.find("Rectangle[")==0){
            if(rect.size()!=4)ABORT("need Rectangle[realMin,realMax,imagMin,imagMax], got: "+Kind);
            lim={tools::string_to_double(rect[0]),tools::string_to_double(rect[1]),
                 tools::string_to_double(rect[2]),tools::string_to_double(rect[3])};
        }
        else if(Kind.find("Range[")==0){
            if(rect.size()<2 or 3<rect.size())ABORT("need Range[realMin,realMax{,nMax}], got: "+Kind);
            lim={tools::string_to_double(rect[0]),tools::string_to_double(rect[1]),-DBL_MAX,DBL_MAX};
        }

        for(int k=_eigenvalues.size()-1;k>=0;k--){
            double eR=_eigenvalues[k].real(),eI=_eigenvalues[k].imag();
            if(lim[0]>eR or eR>lim[1] or lim[2]>eI or eI>lim[3]){
                if(_rightVectors.size()==_eigenvalues.size()){
                    delete _rightVectors[k];
                    _rightVectors.erase(_rightVectors.begin()+k);
                }
                if(_dualVectors.size()==_eigenvalues.size()){
                    delete _dualVectors[k];
                    _dualVectors.erase(_dualVectors.begin()+k);
                }
                if(_leftVectors.size()==_eigenvalues.size()){
                    delete _leftVectors[k];
                    _leftVectors.erase(_leftVectors.begin()+k);
                }
                _eigenvalues.erase(_eigenvalues.begin()+k);
            }
        }
        sort("SmallReal");
    }

    else
        ABORT("selection kind not defined: "+Kind);
    PrintOutput::DEVmessage(Sstr+"selected"+_eigenvalues.size()+"eigenvalues of"+eSize+"computed for problem size"+_ovr->iIndex->sizeStored()+mess);
    LOG_POP();
}

void srtByInt(const std::vector<int> Srt, std::vector<Coefficients*>& C){
    if(!C.size())return;
    if(Srt.size()!=C.size())DEVABORT("wrong sizes");
    std::vector<Coefficients*> cSrt;
    for(size_t k=0;k<Srt.size();k++)cSrt.push_back(C[Srt[k]]);
    std::swap(cSrt,C);
}

void EigenSolverAbstract::sort(string Kind){
    // (really dumb)
    std::vector<int> srt;
    for(size_t k=0;k<_eigenvalues.size();k++)srt.push_back(k);
    tools::sortKey(Kind,_eigenvalues,srt);
    srtByInt(srt,_leftVectors);
    srtByInt(srt,_rightVectors);
    srtByInt(srt,_dualVectors);
}

void EigenSolverAbstract::phaseFix(){
    for(size_t k=0;k<_rightVectors.size();k++){
        // phase at maximum coefficient
        std::complex<double> phaseAtMax=_rightVectors[k]->cMaxNorm();
        phaseAtMax/=std::abs(phaseAtMax);
        _rightVectors[k]->scale(std::conj(phaseAtMax));
        _dualVectors[k]->scale(phaseAtMax);
    }
}

void EigenSolverAbstract::normalize(std::vector<Coefficients*> & Dual, std::vector<Coefficients*> & RVec){
    for(size_t k=0;k<RVec.size();k++){
        // phase at maximum coefficient
        std::complex<double> phaseAtMax=RVec[k]->cMaxNorm();
        phaseAtMax/=std::abs(phaseAtMax);
        // L2-normalize rhs
        double realNrm=std::abs(RVec[k]->idx()->overlap()->matrixElement(*RVec[k],*RVec[k]));
        RVec[k]->scale(std::conj(phaseAtMax)/sqrt(realNrm));
        Dual[k]->scale(1./Dual[k]->innerProduct(RVec[k],true));
    }
}
void EigenSolverAbstract::normalize(){
    normalize(_dualVectors,_rightVectors);
//    for(size_t k=0;k<_rightVectors.size();k++){
//        // phase at maximum coefficient
//        std::complex<double> phaseAtMax=_rightVectors[k]->cMaxNorm();
//        phaseAtMax/=std::abs(phaseAtMax);
//        // L2-normalize rhs
//        double realNrm=std::abs(_ovr->matrixElement(*_rightVectors[k],*_rightVectors[k]));
//        _rightVectors[k]->scale(std::conj(phaseAtMax)/sqrt(realNrm));
//        _dualVectors[k]->scale(1./_dualVectors[k]->innerProduct(_rightVectors[k],true));
//    }
}

void EigenSolverAbstract::orthonormalize(){
    rightVectors();
    dualVectors();
    orthonormalize(_dualVectors,_rightVectors);
}

static size_t verifyLimit=10000000; // at this number it takes a few seconds on my laptop
bool EigenSolverAbstract::verify(double Epsilon) const {
    if(tRecX::off("EigenSolver"))return true;
    if(_epsVerified<=Epsilon)return true;

    // verification is costly for large problems (scales as  size^2*N)
    size_t sample=sqrt(double(std::pow(_rightVectors.size(),2)*_op->idx()->size()/verifyLimit))+1;
    if(ReadInput::main.flag("DEBUGeigenFullVerify",
                            "force complete verification of eigenvectors even when problem is large"))sample=1;

    if(sample>1)PrintOutput::DEVmessage(Sstr+"EigenSolver: large problem - verify sampled at steps of "+sample);

    _epsVerified=0.;
    vector<double> err;
    std::complex<double> errEval;
    for(size_t k=0;k<_rightVectors.size();k+=sample){
        _ovr->apply(1.,*_rightVectors[k],0.,*_ovr->tempLHS());
        _op->apply(1.,*_rightVectors[k],0.,*_op->tempLHS());
        _op->tempLHS()->makeContinuous();
        _op->tempLHS()->axpy(-_eigenvalues[k],_ovr->tempLHS());
        err.push_back(abs(_op->tempLHS()->norm()/(_rightVectors[k]->norm()*max(abs(_eigenvalues[k]),1.))));
        if(_epsVerified<err.back()){
            _epsVerified=err.back();
            errEval=_eigenvalues[k];
        }
    }

    if(_epsVerified>Epsilon){
        PrintOutput::DEVwarning(Sstr+"largest eigenvector error"+_epsVerified+" at eigenvalue"+errEval);
    }

    Eigen::MatrixXcd one=Eigen::MatrixXcd::Identity(std::max(1,int(_dualVectors.size()/sample)+1),std::max(1,int(_rightVectors.size()/sample)+1));
    for(size_t k=0;k<_dualVectors.size();k+=sample)
        for(size_t l=0;l<_rightVectors.size();l+=sample)
            one(k/sample,l/sample)=_dualVectors[k]->innerProduct(_rightVectors[l],true);
    if(not EigenTools::purge(one).isIdentity(Epsilon)){
        one-=Eigen::MatrixXcd::Identity(one.rows(),one.cols());
        PrintOutput::DEVwarning(Str("incorrect duals, max error = ")+one.lpNorm<Eigen::Infinity>());
    }

    if(_epsVerified<Epsilon)_epsVerified=Epsilon*0.999;
    return _epsVerified<Epsilon;
}
