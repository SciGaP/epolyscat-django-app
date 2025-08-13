// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "discretizationSpectral.h"

#include <vector>
#include <complex>

#include "basisAbstract.h"
#include "basisIntegrable.h"

#include "discretizationDerived.h"
#include "indexNew.h"
#include "operatorDiagonal.h"
#include "parallel.h"
#include "operatorFloor.h"

#include "printOutput.h"
#include "readInput.h"
#include "coefficients.h"
#include "operatorAbstract.h"
#include "operatorTree.h"
#include "operatorDefinition.h"
#include "eigenSolver.h"
#ifdef _USE_HACC_
#include "discretizationHaCC.h"
#endif
#include "inverse.h"
#include "tRecXchecks.h"
#include "indexProd.h"
#include "basisSub.h"

#include "parallelOperator.h"

#include "plot.h"
#include "projectSubspace.h"
#include "tools.h"
#include "eigenTools.h"

using namespace std;

DiscretizationSpectral::~DiscretizationSpectral(){
    delete spectralValues;
    spectralValues=0;
}

DiscretizationSpectral::DiscretizationSpectral(const Discretization *D, std::string Criterion)
    : spectralValues(0),_selectionCriterion(Criterion){parent=D;}

DiscretizationSpectral::DiscretizationSpectral(const Discretization *D, const OperatorAbstract *Op,
                                               double Emin, double Emax, bool excludeEnergyRange)
    :DiscretizationSpectral(Op,Emin,Emax,excludeEnergyRange)
{
    parent=D;
    name=D->name+"_"+Op->name+_selectionCriterion;
    name+="_spectral";
}
const std::vector<std::complex<double>>& DiscretizationSpectral::eigenvalues() const{
    return spectralOper()->diagonal();
}

TIMER(dspec,);
TIMER(dspec1,slv);
TIMER(dspec2,slv);
TIMER(dspec3,slv);
TIMER(dspec4,slv);
TIMER(dspec5,slv);
DiscretizationSpectral::DiscretizationSpectral(const OperatorAbstract *Op,
                                               double Emin, double Emax, bool excludeEnergyRange)
    :DiscretizationSpectral(0,"["+tools::str(Emin,3,DBL_MAX/2.)+","+tools::str(Emax,3,DBL_MAX/2.)+"]")
{

    LOG_PUSH("spec1");
    _selectionCriterion="["+tools::str(Emin,3,DBL_MAX/2.)+","+tools::str(Emax,3,DBL_MAX/2.)+"]";
    name=Op->name+"_"+Op->iIndex->hierarchy()+_selectionCriterion+"_spectral";
    idx()=0;

    // solve eigenproblem
    EigenSolver slv(Emin,Emax,true,true,excludeEnergyRange,"Lapack");
    slv.parallel(true);
    bool doFull=true;
    if(doFull)slv.fullVectors();
    slv.compute(Op);

    LOG_POP();
    LOG_PUSH("spec1a");

    if(slv.eigenvalues().size()==0)return;

    // debug
    if(MPIwrapper::Size(MPIwrapper::worldCommunicator())==1 and not slv.verify()){
        UseMatrix mat,ovr;
        Op->matrix(mat);
        mat.print("Op",0);

        Op->iIndex->overlap()->matrix(ovr);
        ovr.print("ovr",0);
        ABORT("Eigensolving failed for "+name);
    }
    LOG_POP();

    LOG_PUSH("spec2");

    if(doFull)
        _project.reset(new ProjectSubspace(slv.rightVectors(),slv.dualVectors()));
    else
        _project.reset(new ProjectSubspace(slv));

    LOG_POP();
    LOG_PUSH("spec3");

    // create index and maps
    std::shared_ptr<OperatorTree> _mapFrom,_mapTo;
    _mapFrom=_project->mapFrom();
    _mapTo=_project->mapTo();
    idx()=const_cast<Index*>(_project->subspaceIndex());

    spectralValues=new OperatorDiagonal("eval("+Op->name+")",idx());
    // ProjectSubspace internally sorts the eigenvectors - adjust
    std::vector<std::complex<double> > eSort;
    for(int k: _project->sorting())eSort.push_back(slv.eigenvalues()[k]);

    spectralOper()->add(eSort,idx());

    // Prevent hierarchies Phi1.Phi1.Eta1.Eigen
    if(idx()->childSize() == 1 && idx()->child(0)->axisName() == idx()->axisName()){
        if(_mapFrom->childSize() != 1) ABORT("Unexpected");
        if(_mapTo->childSize() != 1) ABORT("Unexpected");
        if(_mapFrom->child(0)->iIndex != idx()->child(0)) ABORT("Unexpected");
        if(_mapTo->child(0)->jIndex != idx()->child(0)) ABORT("Unexpected");

        idx() = idx()->child(0);
        _mapFrom = std::shared_ptr<OperatorTree>(_mapFrom->child(0));
        _mapTo =  std::shared_ptr<OperatorTree>(_mapTo->child(0));
        DEVABORT("this will break, likely");

        spectralValues = spectralValues->child(0);

        // TODO: Leaks
        idx()->parentRef() = 0;
        _mapFrom->parentRef() = 0;
        _mapTo->parentRef() = 0;
        spectralValues->parentRef() = 0;
    }

    LOG_POP();
    LOG_PUSH("spec4");

    setToParent(std::shared_ptr<OperatorAbstract>(_mapTo));
    setFromParent(std::shared_ptr<OperatorAbstract>(_mapFrom));
    check(Op);

    Parallel::setSort(this);
    LOG_POP();
}

int DiscretizationSpectral::selectNmax(const Index *Idx) const{
    if(_selectionCriterion=="" or _selectionCriterion[0]=='[')return INT_MAX;

    if(_selectionCriterion.find("Rn<=")!=string::npos){
        // special case for smooth extrapolation
        if(Idx->axisName()!="Rn")ABORT("axis does not match slection "+Idx->axisName());
        double rmax=tools::string_to_double(_selectionCriterion.substr(_selectionCriterion.find("Rn<=")+4));
        int nmax=0;
        for(size_t k=0;k<Idx->childSize() and Idx->child(k)->basis()->integrable()->upBound()*1.0000001<rmax;k++)
            nmax+=Idx->child(k)->basis()->size();
        return nmax;
    }
    else ABORT("undefined eigenvalue selection: "+_selectionCriterion);
    return 0;
}

OperatorTree * DiscretizationSpectral::_mapConstructor(bool ToSpectral, const Discretization* Disc, const Index *SIndex, const Index *EIndex,
                                                       const std::vector<Coefficients*> & Evec)
{
    OperatorTree * mapT;
    if(ToSpectral)mapT=new OperatorTree(Disc->name+" <-- "+Disc->parent->name,SIndex,EIndex);
    else          mapT=new OperatorTree(Disc->parent->name+" <-- "+Disc->name,EIndex,SIndex);

    if(Evec.size()==0)return mapT;

    if(EIndex!=Evec[0]->idx())ABORT("EIndex must match Evec.idx(), is\n"+EIndex->str()+"\n"+Evec[0]->idx()->str());

    // descend to coefficient floor level
    vector<Coefficients*>evec(Evec.size());
    for(unsigned int k=0;k<Evec[0]->childSize();k++){
        for(unsigned int l=0;l<Evec.size();l++)evec[l]=Evec[l]->child(k);
        mapT->childAdd(_mapConstructor(ToSpectral,Disc,SIndex,evec[0]->idx(),evec));
    }

    if(Evec[0]->isLeaf()){
        string hash=name+tools::str(Evec[0]->levelRank());
        vector<complex<double> > mat;
        for(unsigned int k=0;k<Evec.size();k++)
            for(unsigned int l=0;l<Evec[0]->size();l++)
                mat.push_back(Evec[k]->floorData()[l]);

        UseMatrix fmat;
        if(ToSpectral) {
            // to spectral from original
            fmat=UseMatrix::UseMap(mat.data(),Evec[0]->size(),Evec.size()).transpose();
            fmat.purge(1.e-14,1.e-15);
            mapT->floor()=OperatorFloor::factory(vector<const UseMatrix*>(1,&fmat),hash+SIndex->hash()+Evec[0]->idx()->hash());
        }

        else {
            // from spectral to original
            fmat=UseMatrix::UseMap(mat.data(),Evec[0]->size(),Evec.size());
            fmat.purge(1.e-14,1.e-15);
            mapT->floor()=OperatorFloor::factory(vector<const UseMatrix*>(1,&fmat),hash+Evec[0]->idx()->hash()+SIndex->hash());
        }
    }
    return mapT;
}

void DiscretizationSpectral::check(const OperatorAbstract* Op) const {
    if(tRecX::off("spectralMap"))return;
    if(idx()==0 or idx()->size()==0){
        PrintOutput::DEVwarning("cannot check empty discretization: "+name);
        return;
    }

    Str mess(""," ");
    Coefficients cspec(idx());
    cspec.setToRandom();
    cspec*=1.;
    Coefficients cinit(cspec);
    Coefficients cparent(Op->iIndex);
    mapToParent()->apply(1.,cspec,0.,cparent);
    mapFromParent()->apply(-1.,cparent,1.,cinit);
    if(not cinit.isZero(1.e-9)){
        mess=mess+"\nnot identity on spectral space:\n"+cinit.str(2);
    }
    Coefficients caux(cparent);
    mapFromParent()->apply(1.,cparent,0.,cinit);
    mapToParent()->apply(-1.,cinit,1.,caux);
    if(not caux.isZero(1.e-9)){
        mess=mess+"\nnot projector on parent space:\n"+caux.str();
    }

    double errMax=0.;
    if(Op != 0){
        Coefficients c1(Op->jIndex);
        Coefficients c1Proj(Op->jIndex);
        Coefficients c2(Op->iIndex);
        Coefficients cSpec(idx());
        Coefficients cSpec2(idx());
        Coefficients c2Check(Op->iIndex);

        c1.setToRandom();
        for(size_t k=0;k<c1.size();k++)c1.data()[k]=k;
        c1.makeContinuous();

        mapFromParent()->apply(1., c1, 0., cSpec);
        mapToParent()->apply(1., cSpec, 0., c1Proj);

        Op->apply(1., c1Proj, 0., c1);
        Op->iIndex->inverseOverlap()->apply(1., c1, 0., c2);
        c2.makeContinuous();

        spectralOper()->updateFunction(0., OperatorDiagonal::identityFunction);
        spectralOper()->apply(1., cSpec, 0., cSpec2);
        mapToParent()->apply(1., cSpec2, 0., c2Check);

        // get error, relative to element sizes
        double eps=1e-7;
#ifdef _DEVELOP_
        eps=1e-10;
#endif
        if(not c2.cwiseRelativeError(c2Check).isZero(eps)){
            mess=mess+"HP=U^\\dagger d U not satisfied, error=";
            mess=mess+c2.norm();
        }
        errMax=std::max(c2.norm(),errMax);
    }

    if(mess=="")PrintOutput::DEVmessage(Str("OK spectral maps for")+name+"(size"+cinit.size()+")");
    else        PrintOutput::warning(mess+"for"+name+"(size"+cinit.size()+")",1,0,
                                     " This is a test on the spectral decomposition"
                                     "\n Numerical errors are inevitable, but the level appears high and may compromize results"
                                     "\n Possible fix: lower order in discretization, higher cutoff energy"
                                     );
    if(errMax>1.e-6){
        Sstr+"Op\n"+Op->str()+Sendl;
        Sstr+"Op.idx\n"+Op->iIndex->str()+Sendl;
        Sstr+"eigenvalues\n"+spectralOper()->str()+Sendl;
        DEVABORT("severe spectral error");
    }
}
























