// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "initialState.h"

#include "mpiWrapper.h"
#include "readInput.h"
#include "wavefunction.h"
//#include "operator.h"
#include "discretizationConstrained.h"
#include "derivativeFlat.h"

#include "printOutput.h"
#include "eigenSolver.h"
#include "eigenSolverAbstract.h"
#include "plot.h"
#include "parallelOperator.h"
#include "basisOrbital.h"
#include "vectorValuedFunction.h"
#include "productFunction.h"
#include "basisOrbitalNumerical.h"
#include "basisEigen.h"
#include "tools.h"
#include "eigenSolverNonLin.h"
#include "checkpoint.h"

using namespace std;

InitialState::InitialState(){}

std::string readPrint();

std::string InitialState::readKind(ReadInput & Inp){
    std::string kind;
    Inp.texdocuCategoryAdd("Initial","kind,state","","60,80,220");
    Inp.read("Initial","kind",kind,"atBegin",
             string("initial state: atBegin...eigenstate of time-evolution at t=tBegin")+
             string(", mln...for polar select by quantum numbers, Hinitial...eigenstate of specified operator")+
             string(", inFlux...zero, incoming flux, manyBody...many-body calc")+
             string(", Orbital...orbital in basis")+
             string(", ProductFunction...replaces Function")+
             string(", External...numerical external orbital")+
             string(", Eigenstate...[oper:n:refIdx] n'th eigenstate of operator wrt refIdx")+
             string(", EigenSubblock...[Hinitial(i0,i1,..):n] n'th eigenstate of (i0,i1,..) diagonal subblock of Hinitial")+
             string(", GP...ground state of Gross Pitaeskii")+
             string(", EE...mean field electron-electron potential")
             )
            .texdocu(R"tex(
                     How to construct $\Psi_0=\Psi(t_0)$ for time-propagation starting at $t_0$:
                     \begin{itemize}
                     \item \lcode{atBegin}: lowest eigenstates of $H(t_0)$ [Operator:\nameref{docu:Operator:hamiltonian}]
                     \item \lcode{EigenSubblock} $n$'th eigenvector $\phi_n$ of a diagonal subblock of the full
                     Operator:\nameref{docu:Operator:initial}. $\Psi_0=(0,\ldots,0,\phi_n,\ldots)$,
                     where $\phi_n$ is in the position of the subblock. The subbblock is selected by the multiindex
                     \lcode{i0,i1,..}\\
                     Simple example: for the \lcode{i0,i1}=$m,l$ angular quantum numbers,
                     this selects the $n$'th excited state in the m,l symmetry.
                     \\
                     In general, \lcode{Hinitial} does not need to be block-diagonal.
                     \item \lcode{Eigenstate} For use in hybrid bases: $\Psi_0$ is $n$'th eigenvector \wrt on \nameref{docu:Axis:subset}.
                     \lcode{operator} can be any legitimate operator string for \lcode{refIdx}. Used, e.g., in \nameref{docu:tutorial:220}
                     \item \lcode{mln} [OBSOLESCENT]
                     \item \lcode{inFlux} [DEVELOPMENT]
                     \end{itemize}
                     )tex");
    return kind;
}

// transform Vecs such that they solve the eigenproblem on the unscaled region
void diagonalizeUnscaled(const OperatorTree & Ham,const OperatorTree & Ovr,
                         const vector<Coefficients*> & In,vector<Coefficients> & Out){
    Eigen::MatrixXcd h(In.size(),In.size());
    Eigen::MatrixXcd o(In.size(),In.size());
    Coefficients hamV(Ham.iIndex),ovrV(Ovr.iIndex),tmp(Ovr.jIndex);
    for(int k=0;k<In.size();k++){
        tmp=*In[k]; // Ham only matches first block of Vecs
        Ham.apply(1.,tmp,0.,hamV);
        Ovr.apply(1.,tmp,0.,ovrV);
        for (int l=0;l<=k;l++){
            h(l,k)=In[l]->innerProductUnscaled(&hamV);
            o(l,k)=In[l]->innerProductUnscaled(&ovrV);
            h(k,l)=conj(h(l,k));
            o(k,l)=conj(o(l,k));
        }
    }

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> slv;
    slv.compute(h,o);
    Out.clear();
    for(int k=0;k<In.size();k++){
        Out.push_back(Coefficients(Ham.iIndex));
        for(int l=0;l<In.size();l++)
            Out[k].axpy(slv.eigenvectors()(l,k),In[l]);
    }
}

TIMER(get,)
Wavefunction InitialState::get(const Discretization* D, std::string initialKind, int initialN,
                               double tBeg,
                               const OperatorTree &InitialOper,
                               DerivativeFlat *derNew){




    std::string initialStateMessage="initial state: "+initialKind;
    Coefficients cIni(D->idx());
    Wavefunction wf(tBeg,&cIni);
    wf.coefs->setToZero();

    // use checkpointed (if available)
    Checkpoint chPt(ReadInput::main.output());
    if(chPt()){
        wf.time=chPt.time();
        chPt.coefs(wf.coefs);
        PrintOutput::message(Sstr+"state from checkpoint at time"+wf.time);
        return wf;
    }

    START(get);
    vector<complex<double> > Eval;
    vector<Coefficients*>Evec;
    unsigned int nini=max(initialN,1)-1;

    if(initialKind=="Hinitial"){
        // compute as defined separately on input
        // start from H0 ground state
        if(not DiscretizationConstrained::inputs(ReadInput::main)){
            EigenSolver slv(-DBL_MAX,DBL_MAX,true,false,false);
            slv.compute(&InitialOper);
            slv.select("SmallReal["+tools::str(nini+1)+"]");
            std::vector<Coefficients*> vecs = slv.rightVectors();
            vecs.back()->scale(1./sqrt(vecs.back()->idx()->overlap()->matrixElement(*vecs.back(), *vecs.back())));
            wf=Wavefunction(tBeg,vecs.back());
        }
        else {
            //NOTE: should be converted to OperatorTree
            DiscretizationConstrained constrD(D,ReadInput::main);
            DEVABORT("Constraint input has not been used in a long time, fix and test");
//            Operator constrH("constrainedH",InitialOper.def(),&constrD,&constrD); // hamiltonian for initial state
//            constrH.eigen(*constrH.iIndex->localOverlap(),Eval,Evec,nini+1);
//            constrD.mapToParent()->apply(1.,*Evec[nini],0.,cIni);
//            cIni.scale(1./sqrt(cIni.idx()->localOverlap()->matrixElement(cIni,cIni)));
//            wf=Wavefunction(tBeg,&cIni);
        }
        initialStateMessage="initial state: eigenvector "+tools::str(nini+1)+" of operator "+InitialOper.name;

        // make sure the starting wave function is properly projected
        CoefficientsGlobal* globView=CoefficientsGlobal::view(wf.coefs);
        complex<double> nrm0=wf.coefs->idx()->localOverlap()->matrixElement(*globView,*globView);
        derNew->project(*globView);
        complex<double> nrm1=wf.coefs->idx()->localOverlap()->matrixElement(*wf.coefs,*wf.coefs);
        if(abs(nrm1-nrm0)>abs(nrm0)*1.e-8){
            PrintOutput::warning("projection changed initial norm of "+tools::str(nrm0)+" by "+tools::str(nrm1-nrm0));
            PrintOutput::warning("projection may be too ambitious, increase TimePropation:cutEnergy");
        }
    }
    else if(initialKind=="atBegin"){
        derNew->update(tBeg);
        derNew->eigen(Eval,Evec,nini+1);
        Eval.clear();

        if(MPIwrapper::isMaster())wf=Wavefunction(tBeg,Evec[nini]);
        initialStateMessage="initial state: eigenvector "+tools::str(nini+1)+" of time-evolution at t="+tools::str(tBeg);
    }

    else if(initialKind=="manyBody"){
        // the lowest single particle eigenvectors are the components of the many-body initial state
        const OperatorTree* op=&InitialOper;
        if(op->name!="H0")
            ABORT("for Initial:kind=many define Operator: initial=Hamiltonian, is "+op->name);
        if(op->iIndex->axisName()!="Vec")
            ABORT("not a many-body discretization, hierarchy="+op->iIndex->hierarchy());
        while(op->iIndex->axisName()=="Vec")op=op->child(0);
        EigenSolver iniSlv(-DBL_MAX,1.);
        iniSlv.compute(op);
        iniSlv.sort("SmallReal");
        PrintOutput::title("SINGLE PARTICLE EIGENVALUES");
        PrintOutput::newRow();
        PrintOutput::rowItem("real(E)");
        PrintOutput::rowItem("imag(E)");
        if(iniSlv.eigenvalues().size()<wf.coefs->childSize())
            ABORT(Str("too few eigenvalues < 0:")+iniSlv.eigenvalues().size()+"with values"+iniSlv.eigenvalues());
        for(int k=0;k<min((int)iniSlv.eigenvalues().size(),(int)wf.coefs->childSize()+3);k++){
            PrintOutput::newRow();
            PrintOutput::rowItem(real(iniSlv.eigenvalues()[k]));
            PrintOutput::rowItem(imag(iniSlv.eigenvalues()[k]));
            if(k<wf.coefs->childSize())PrintOutput::rowItem("*");

        }

        // get eigenvectors wrt unscaled region
        vector<Coefficients*>inVec;
        vector<Coefficients> outVec;
        for(int k=0;k<wf.coefs->childSize();k++)inVec.push_back(iniSlv.rightVectors()[k]);
        diagonalizeUnscaled(op,op->iIndex->overlap(),inVec,outVec);
        for(int k=0;k<wf.coefs->childSize();k++)*wf.coefs->child(k)=outVec[k];

        // normalize total to 1
        double nrm=real(wf.coefs->idx()->overlap()->matrixElement(*wf.coefs,*wf.coefs));
        wf.coefs->scale(1./sqrt(nrm));
        PrintOutput::paragraph();
    }
    else if(initialKind.find("Function[")==0){
        ABORT("OBSOLETE "+initialKind+" use ProductFunction[...] instead");
    }
    else if(initialKind=="EE"){
        /*
        const OperatorTree* op=&InitialOper;
        Coefficients psiZero(op->jIndex);
        psiZero.setToZero();


        derNew->update(tBeg, &psiZero);
        derNew->eigen(Eval,Evec,nini+1);
        Eval.clear();
        if(MPIwrapper::isMaster())wf=Wavefunction(tBeg,Evec[nini]);
        initialStateMessage="initial state: eigenvector "+tools::str(nini+1)+" of time-evolution at t="+tools::str(tBeg);

*/

        const OperatorTree* op=&InitialOper;
        Coefficients psiZero(op->jIndex);
        psiZero.setToZero();
        derNew->update(tBeg, &psiZero);
        EigenSolverNonLin slv(-DBL_MAX,DBL_MAX,INT_MAX,true,false,false,"Lapack");
        slv.compute(op);
        slv.sort("SmallReal");
        wf=Wavefunction(tBeg,slv.rightVectors()[0]);


    }
    else if(initialKind=="GP"){
        derNew->update(tBeg);
        const OperatorTree* op=&InitialOper;
        EigenSolverNonLin slv(-DBL_MAX,DBL_MAX,INT_MAX,true,false,false,"Lapack");
        slv.compute(op);
        slv.sort("SmallReal");
        wf=Wavefunction(tBeg,slv.rightVectors()[0]);
    }

    if(wf.coefs->isZero())productFunction(initialKind,*wf.coefs);
    if(wf.coefs->isZero())indexOrbital(initialKind,*wf.coefs);
    if(wf.coefs->isZero())externalOrbital(initialKind,*wf.coefs);
    if(wf.coefs->isZero())subblockEigen(initialKind,&InitialOper,*wf.coefs);

    if(MPIwrapper::isMaster() and initialKind!="ZERO" and wf.coefs->isZero())
        ABORT("not implemented: Initial: kind="+initialKind);

    PrintOutput::message(initialStateMessage);
    STOP(get);
    return wf;
}

void InitialState::productFunction(string InitialKind, Coefficients& C){
    C.setToZero();
    if(InitialKind.find("ProductFunction[")!=0)return;

    vector<string> part=tools::splitString(tools::stringInBetween(InitialKind,"[","]"),':');
    vector<string> prod=tools::splitString(part[1],',',"<{[(",">}])");
    std::shared_ptr<VectorValuedFunction>func(new ProductFunction(part[0],prod));
    VectorValuedFunction::add("InitialState",func);
    BasisOrbitalNumerical coefs("InitialState:main:0:1");
    C=*coefs.orbital(0);
}

void InitialState::indexOrbital(string InitialKind, Coefficients& C){
    C.setToZero();
    if(InitialKind.find("Orbital[")!=0)return;

    Coefficients* c=&C;
    while(c!=0 and c->idx()->axisName()!="Orbital")c=c->nodeNext();
    if(c==0)ABORT("Cannot use Initial:Orbital, no Orbital axis in hierarchy "+C.idx()->hierarchy());

    const BasisOrbital* b=dynamic_cast<const BasisOrbital*>(c->idx()->basis());
    if(b==0)ABORT("axis Orbital does not have Orbital basis: "+c->idx()->strNode());
    int iOrb=tools::string_to_int(tools::stringInBetween(InitialKind,"Orbital[","]",true));
    c->orderedData()[iOrb]=1.;
    complex<double> nrm=C.idx()->overlap()->matrixElement(C,C,false);
    c->scale(1./sqrt(nrm));

    Plot plt(b->orbital(iOrb)->idx(),ReadInput::main);
    if(not plt.isEmpty()){
        std::string pltFile=ReadInput::main.output()+"InitialState";
        plt.plot(*b->orbital(iOrb),pltFile);
        PrintOutput::message(Sstr+"Initial state on"+pltFile+"(orbital"+iOrb+"of"+b->name()+")");
    }

}

void InitialState::externalOrbital(string InitialKind, Coefficients& C){
    C.setToZero();
    if(InitialKind.find("External[")!=0 and InitialKind.find("Eigenstate[")!=0)return;
    vector<string> part=tools::splitString(tools::stringInBetween(InitialKind,"[","]"),':');
    int iStat=0;
    if(part.size()>1)iStat=tools::string_to_int(part[1]);
    std::string refIdx="main";
    int instance=0;
    if(part.size()>2){
        refIdx=part[2].substr(0,part[2].find("{"));
        if(part[2].find("{")!=string::npos){
            if(part[2].find("{all}")!=string::npos)
                instance=-1;
            else
                instance=tools::string_to_int(tools::stringInBetween(part[2],"{","}"));
        }
    }
    std::unique_ptr<BasisOrbital> orbs;
    if(InitialKind.find("External[")==0){
        BasisOrbitalNumerical * b=new BasisOrbitalNumerical(part[0]+":"+refIdx+":"+tools::str(iStat)+":1");
        b->generateOrbitals();
        orbs.reset(b);
        b->print("Initial state");
    }
    else if(InitialKind.find("Eigenstate[")==0){
        orbs.reset(new BasisEigen(part[0]+":"+refIdx,1,iStat));
        PrintOutput::message(Str("selected orbital"," ")+orbs->str());
    }

    //  place into desired branch
    Coefficients* subC=&C;
    if(refIdx!="main"){
        for(;subC!=0 and subC->idx()->axisSubset()!=refIdx;subC=subC->nodeNext());
        if(subC==0)ABORT("no reference index \""+refIdx+"\" in main\n"+C.idx()->str());
    }

    double nrm= subC->parent()==0 or instance!=-1 ? 1. : 1./sqrt(double(subC->parent()->childSize()));
    for(int k=0;subC!=0;k++,subC=subC->nodeRight()){
        if(subC==0)ABORT(Str("no branch ","")+instance+" with reference \""+refIdx+"\" in\n"+C.idx()->str());
        if(k==instance or instance==-1)*subC=*orbs->orbital(0);
    }
    C*=nrm;
}

void InitialState::subblockEigen(string InitialKind, const OperatorTree* Op, Coefficients& C){
    C.setToZero();
    if(InitialKind.find("EigenSubblock[")!=0)return;

    std::vector<std::string> part=tools::splitString(tools::stringInBetween(InitialKind,"[","]"),':');
    if(part[0].find("Hinitial(")!=0)
        ABORT("limited to Hinitial for now, specify as \"Hinitial(2,3)\" for the diagonal block 2 and 3 on first two hierarchy levels, found: "+part[0]);
    std::vector<std::string> sBlock=tools::splitString(tools::stringInBetween(part[0],"(",")"),',');
    if(sBlock[0].find("Hinitial")==0)
        ABORT("specify block as \"Hinitial(2,3)\" for the diagonal block 2 and 3 on first two hierarchy levels, found: "+part[0]);
    std::vector<unsigned int>block;
    for(std::string s: sBlock)block.push_back(tools::string_to_int(s));

    int iStat= part.size()<2 ? 0 : tools::string_to_int(part[1]);

    // descend to diagonal block
    const OperatorTree* op=Op;
    while(op!=0 and op->iIndex->index()!=block and op->jIndex->index()!=block)op=op->nodeNext(Op);
    if(op==0)ABORT(Sstr+"diagonal block"+block+"not found in"+Op->name+"\nDefinition\n"+Op->def());

    EigenSolver slv(-100.,100.,iStat+1,true,true,false,"Arpack");
    slv.compute(op);
    slv.sort("SmallReal");
    *C.nodeAt(block)=*slv.rightVectors()[iStat];
}



