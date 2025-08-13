// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "qtEigenSparse.h"
#include "operatorTree.h"
#include "tools.h"
#include "timeCritical.h"
#include "mpiWrapper.h"
#include "tRecX_cache.h"

#include "discretizationHybrid.h"

#include "printOutput.h"
#include "index.h"
#include "coefficients.h"
#include "operatorFloor.h"
#include "operatorFloorNonLin.h"
#include "operatorData.h"
#include "parameters.h"
#include "operatorTensor.h"
#include "parallel.h"
#include "str.h"
#include "basisNdim.h"
#include "diagnose.h"

#include "operatorZG.h"
#include "operatorZD.h"
#include "operatorFloorEE.h"
#include "operatorZDxZG.h"
#include "operatorZGxZD.h"
#include "operatorZGxZG.h"
#include "operatorZGSblock.h"
#include "constrainedView.h"

#include "operatordDiagsPermuted.h"
#include "operatorPermute.h"

#include "asciiFile.h"
#include "basisMat1D.h"
#include "basisGrid.h"

#include "parallelOperator.h"

#include "log.h"
#include "algebra.h"
#include "basisMatDependent.h"
#include "basisOrbital.h"
#include "basisOrbitalNumerical.h"
#include "basisExpIm.h"
#include "basisChannel.h"

#include "operatorFloorXC.h"
#include "operatorHartree.h"

#include "parallelOperator.h"
#include "plot.h"

using namespace std;

bool OperatorTree::debug=false;

static double epsPurge=1.e-12;

static std::map<std::string,size_t> memoryUsage;

static std::map<const Index*, size_t> indexSections(const Index* Idx, size_t MaxDepth){
    std::map<const Index*, size_t> res;
    int iSec=0;
    for(const Index* ix=Idx;ix and ix->depth()<=MaxDepth;ix=ix->nodeNext()){
        if(ix->depth()==MaxDepth or ix->hasFloor())res[ix]=iSec++;
    }
    return res;
}

Eigen::MatrixXd OperatorTree::matrixNorms(size_t MaxDepth) const {
    auto iSecs=indexSections(idx(),MaxDepth);
    auto jSecs=indexSections(jdx(),MaxDepth);
    Eigen::MatrixXd res(Eigen::MatrixXd::Zero(iSecs.size(),jSecs.size()));
    for(const OperatorTree* op=this;op;op=op->nodeNext()){
        if(iSecs.count(op->idx()) and iSecs.count(op->jdx()))
            res(iSecs[op->idx()],jSecs[op->jdx()])=op->matrix().lpNorm<Eigen::Infinity>();
    }
    return res;
}

OperatorTree::~OperatorTree(){
    ParallelOperator::unsetHost(this);
    if(not _view)delete oFloor;
    oFloor=0;
}

UseMatrix* OperatorTree::colMat;
std::vector<std::complex<double> > OperatorTree::colStor;

OperatorTree::OperatorTree(const std::string Name, const std::string Definition, const Index* IIndex, const Index* JIndex, OperatorFloor *OFloor)
    :OperatorAbstract(Name,IIndex,JIndex),oFloor(OFloor),_view(false){definition=Definition;}

TIMER(operFloor,opTree)
TIMER(operPar,opTree)
TIMER(floor2,opTree)
TIMER(fuse,opTree)
TIMERRECURSIVE(opAll,opTree)
TIMER(iSub,opTree)
TIMER(para,opTree)
TIMER(purge,opTree)
TIMER(dist,opTree)
TIMER(parOptree,)
TIMER(fuseOptree,)
TIMER(purgeOptree,)
TIMER(distriOptree,)

static std::map<std::complex<double>*,std::string> functionStrings;
bool OperatorTree::fromCache(std::string Name){
    std::string src=tRecX_cache::source(Name)+"_"+iIndex->hierarchy()+"_"+jIndex->hierarchy();
    if(not folder::exists(src))return false;

    std::ifstream s(src);
    OperatorTree cached(s,name,iIndex,jIndex);

    std::string nameSource=tRecX_cache::source(Name);
    if(definition!=cached.definition)
        ABORT("cannot use cached "+nameSource+", definitions differ: "+definition+ " vs. "+cached.definition);
    if(not cached.iIndex->treeEquivalent(iIndex) or
            not cached.jIndex->treeEquivalent(jIndex))
        ABORT("cannot use cached "+nameSource+" --- indices differ");

    // move all children to present
    cached.childrenMove(*this);
    ParallelOperator par(this);
    par.setDistribution();
    PrintOutput::message("retrieved operator from cache "+nameSource);

    return true;
}


static bool Basis_hasOverlap(const BasisAbstract* IBas,const BasisAbstract* JBas){
    if(IBas->integrable() and JBas->integrable())
        return IBas->integrable()->lowBound()<JBas->integrable()->upBound()
                and JBas->integrable()->lowBound()<IBas->integrable()->upBound();
    if(IBas->grid() and JBas->grid())
        return IBas->grid()->mesh().front()<JBas->grid()->mesh().back()
                and JBas->grid()->mesh().front()<IBas->grid()->mesh().back();
    DEVABORT("Illegal bases with FE: "+IBas->str()+JBas->str());
}


static Eigen::MatrixXcd overFE(const Index* IdxFE,const Index* JdxFE){
    std::vector<const BasisAbstract*> iBas,jBas;
    Eigen::MatrixXcd over(IdxFE->childSize(),JdxFE->childSize());
    for(size_t k=0;k<IdxFE->childSize();k++)
        for(size_t l=0;l<JdxFE->childSize();l++)
            over(k,l)=double(Basis_hasOverlap(IdxFE->child(k)->axisIndex(IdxFE->axisName())->basis(),
                                              JdxFE->child(l)->axisIndex(JdxFE->axisName())->basis()));
    return over;
}

bool isMultiDim(const BasisAbstract* Bas){
#ifdef _USE_HACC_
    return(Bas->ci() or Bas->orbital() or Bas->ndim());
#else
    return(Bas->orbital() or Bas->ndim());
#endif
}

void OperatorTree::addMulti(const string Name, const OperatorDefinition & Definition, const Index *ISub, const Index *JSub,
                            std::complex<double> Multiplier, std::vector<std::complex<double> *> TFac){

    // split into terms, unbracket, expand, such that each term[n] starts as factor<def>... or factor[[defSpecial]]...
    string icoor=ISub->basis()->ndim()?ISub->basis()->ndim()->quadCoor():ISub->coordinates();
    string jcoor=JSub->basis()->ndim()?JSub->basis()->ndim()->quadCoor():JSub->coordinates();
    std::vector<OperatorDefinition> terms;

    if(isMultiDim(ISub->basis()) or isMultiDim(JSub->basis()))
        terms=Definition.terms();
    else
        terms=Definition.singleTerms(icoor,jcoor);

    // if coordinates indicate sub-region, drop terms that are not defined on the sub-region
    if(name!="Commutator"
            and not ISub->basis()->ndim() and not JSub->basis()->ndim()
            and not ISub->basis()->orbital() and not JSub->basis()->orbital()
            and not ISub->basis()->hybrid() and not JSub->basis()->hybrid()
            ){
        OperatorDefinition::dropTerms(terms,ISub->hierarchy());
    }

    // with terms sorted out, construct operator term-wise
    for(auto t: terms){
        if(t.find("undefined")!=std::string::npos)DEVABORT(Sstr+Definition+terms+icoor+jcoor);
        childAdd(new OperatorTree(Name,t,ISub,JSub,Multiplier,TFac));
        // remove if leaf w/o floor (would be purged anyway)
        if(childBack()->childSize()==0 and childBack()->floor()==0)childErase(childSize()-1);
    }
}

OperatorTree::OperatorTree(const std::string Name, const std::string & Definition, const Index* IIndex, const Index* JIndex)
    :OperatorAbstract(Name,IIndex,JIndex),oFloor(0),_view(false)
{
    LOG_PUSH(Name);
    definition=Definition;
    Parameters::updateToOne();
    Timer::generalMonitor.monitor(0,"Operator: "+name);
    addMulti(Name,Definition,IIndex,JIndex,1.,{});

    postProcess();
    Parameters::restoreToTime();
    LOG_POP();
    if(tRecX_cache::source(Name)!=""){
        std::ofstream s(tRecX_cache::source(Name)+"_"+iIndex->hierarchy()+"_"+jIndex->hierarchy(),std::ios::binary);
        write(s);
        PrintOutput::message("Operator cached to "+tRecX_cache::source(Name));
    }
}

OperatorTree::OperatorTree(const string Name, const OperatorDefinition & Definition, const Index *IIndex, const Index *JIndex,
                           std::complex<double> Multiplier, std::vector<std::complex<double> *> TFac)
    :OperatorAbstract(Name,IIndex,JIndex),oFloor(0),_view(false)
    // Recursive construction of OperatorTree:
    //  | (special case: OperatorTree node that involves BasisOrbital
    //  | determine coordinates of present block
    //
    //  | if time-dependent termFactor: add to list
    //  | else: multiply onto Multiplier
    //
    //  | if floor on both indices:
    //      generate floor block
    //  | else if skipAxis:
    //  |   treat as FE-like axis, i.e. check for local and non-local, set tensor factor as Id or AllOnes
    //  |   descend to next left and/or right Index level
    //  | else: determine tensor factor matrix mult:
    //      | standard operators, i.e not containing [[...]]:
    //          | case: mult=<idx().basis()|singleFactor|jdx().basis()>
    //          | case: resolve dependence for factors <{}> (non-tensor product multiplicative factors)
    //      | [[...]], non-standard operators:
    //          |       mult=AllOnes
    //  | sub-Definition for next recursion level:
    //      | if has leading <singleFactor>: replace with remainder
    //      | else:                          leave unchanged
    //  | decide whether or not to descend to sub-Indices
    //  | loop through all sub-Indices:
    //  |     attach OperatorTree children for sub-Defintion with factors mult(isub,jsub)
    //
    // | on entry level: post process
{

    if(idx()->basis()->hybrid()!=jdx()->basis()->hybrid())
        DEVABORT("hybrid must be same on both sides, got: "+idx()->strNode()+Name+jdx()->strNode());

    definition=Definition;

    if(IIndex->parent()==0 and JIndex->parent()==0 and fromCache(Name))return;

    if(OperatorAbstract::useTensor and OperatorAbstract::fuseOp)
        PrintOutput::DEVwarning("simultaneous use of fuse() and true tensor may lead to incorrect results: "+definition,1);

    if(abs(Multiplier)<epsPurge)
        PrintOutput::warning("found near-zero multiplier in setup of operator "+definition);

    // set flag that index build may be ongoing
    Index::build=true;

    STARTDEBUG(opAll);

    OperatorDefinition(Definition).syntax();
    if(IIndex->axisName()=="Orbital" or JIndex->axisName()=="Orbital"){
        if(Definition.find("<<Overlap>>")>Definition.find("<")
                or (IIndex->axisSubset()=="Complement")==(JIndex->axisSubset()=="Complement")){
            BasisOrbital::attachOperator(this,Name,Definition,IIndex,JIndex,Multiplier,TFac);
        }
    }
    else {
        string icoor=iIndex->basis()->ndim()?iIndex->basis()->ndim()->quadCoor():iIndex->coordinates();
        string jcoor=jIndex->basis()->ndim()?jIndex->basis()->ndim()->quadCoor():jIndex->coordinates();
        std::vector<OperatorDefinition> terms=Definition.terms();
        if(terms.size()!=1)DEVABORT(Sstr+"here exactly one term is needed, got"+terms);

        // separate factor and operator
        OperatorDefinition termOper;
        string termFactor=OperatorDefinition(terms[0]).parameter(termOper);

        // add to list of function multipliers or multiply
        if(Parameters::isFunction(termFactor)){
            // use variable factors always with + sign
            termFactor=tools::cropString(termFactor);
            if(termFactor[0]=='-'){
                Multiplier*=-1.;
                termFactor[0]='+';
            } else if(termFactor[0]!='+')
                termFactor="+"+termFactor;
            TFac.push_back(Parameters::pointer(termFactor));
            functionStrings[TFac.back()]=termFactor;
        }
        else
            Multiplier*=*Parameters::pointer(termFactor);

        if((iIndex->hasFloor() and jIndex->hasFloor())) {

            // remove auxiliary setup info (like $BAND.0.0 etc.)
            if(termOper.find("$") != string::npos)
                termOper=termOper.substr(0,termOper.find("$"));

            if(termOper.find("|")!=string::npos)ABORT("inhomogenous term not implemented");


            //PARALL - basics
            // parallelization of operator setup works as follows:
            // the full tree down to floor is built on all nodes
            // floors are only constructed on the floor's "host"
            // the host is determined below: for Hartree and XC, there is a postprocessing
            // step that envolves each floor on a given patch at a time
            // the factory, actually, returns only an empty Hartree or XC place holder
            // these should be on ALL hosts and only in postProcess they will be distributed in memory

            OperatorFloor * of;
            if(definition.find("<Hartree")!=std::string::npos or definition.find("<XC")!=std::string::npos){
                    // needs post-process, set for all (note: this may be moved nt Parallel eventually)
                    of=OperatorFloor::factory(termOper,iIndex,jIndex,Multiplier);
            }
            else if(definition.isStandard(idx(),jdx()) and
                    iIndex->basis()->ndim()==0 and
                    jIndex->basis()->ndim()==0){
                    // set up only on owner
                    of=Parallel::operatorFloor(iIndex, jIndex,
                                               [&](){ return OperatorFloor::factory(termOper,iIndex,jIndex,Multiplier); }
                    );
            }
            else {
                // distribute operator floor and set dummy where empty
                of=Parallel::operatorFloor(iIndex, jIndex,
                                           [&](){
                    auto ret = OperatorFloor::specialFactory(name,termOper,iIndex,jIndex,Multiplier);
                    if(ret == nullptr) ret = new OperatorDUM(0.);
                    return ret;
                }
                );
            }
            if (of==0)DEVABORT("Parallel::operatorFloor returned nullptr");

            // attach time-dependent factor to floor
            if(TFac.size()>0)of->setFactor(TFac[0]);
            if(TFac.size()>1)DEVABORT("multiple function factors no implemented yet");

            definition = termFactor+termOper;
            oFloor = of;
        }
        else{
            // not at floor yet, descend to either side
            const Index *iSub0=iIndex;
            const Index *jSub0=jIndex;
            OperatorDefinition opSub=termOper;
            Eigen::MatrixXcd mult;

            // treat various types of axes as if finite elements level
            if(OperatorDefinition::skipAxis(iIndex) and OperatorDefinition::skipAxis(jIndex)){
                // set up coupling between FE-like sections
                if(iIndex->axisName()!=jIndex->axisName()  and iIndex->axisName().substr(0,4) != "spec")ABORT("situation not covered");
                if(termOper.find("NONORTH")!=std::string::npos)DEVABORT(termOper+": NONORTH no longer supported, use allOnes instead");
                if(not termOper.isLocal(iIndex,jIndex)){
                    mult=Eigen::MatrixXcd::Constant(iIndex->childSize(), jIndex->childSize(), 1.);
                } else {
                    if(iIndex->axisName()!=jIndex->axisName())
                        DEVABORT("FE-axes do not match: "+iIndex->axisName()+jIndex->axisName());
                    if(iIndex->depthOfDuplicate()!=Index::npos and jIndex->depthOfDuplicate()!=Index::npos)
                        // FE axis, check for overlapping elements
                        mult=overFE(iIndex,jIndex);
                    else
                        // FE with equal depth and local operator: different elements do not overlap
                        mult=Eigen::MatrixXcd::Identity(iIndex->childSize(),jIndex->childSize());
                }
                iSub0=iSub0->child(0);
                jSub0=jSub0->child(0);
            }
            else if(OperatorDefinition::skipAxis(iIndex)){
                //                    mult=UseMatrix::Constant(iIndex->childSize(),1,1.);
                mult=Eigen::MatrixXcd::Constant(iIndex->childSize(),1,1.);
                iSub0=iSub0->child(0);
                if(iIndex->axisName()==jIndex->axisName() and iIndex->depthOfDuplicate()!=Index::npos){
                    for(size_t k=0;k<iIndex->childSize();k++)mult(k,0)=
                            double(Basis_hasOverlap(iIndex->child(k)->axisIndex(iIndex->axisName())->basis()->integrable(),
                                                    jIndex->basis()->integrable()));
                }
            }
            else if(OperatorDefinition::skipAxis(jIndex)){
                //                    mult=UseMatrix::Constant(1,jIndex->childSize(),1.);
                mult=Eigen::MatrixXcd::Constant(1,jIndex->childSize(),1.);
                jSub0=jSub0->child(0);
                if(iIndex->axisName()==jIndex->axisName() and jIndex->depthOfDuplicate()!=Index::npos){
                    for(size_t k=0;k<jIndex->childSize();k++)mult(0,k)=
                            double(Basis_hasOverlap(jIndex->child(k)->axisIndex(jIndex->axisName())->basis()->integrable(),
                                                    iIndex->basis()->integrable()));
                }
            }
            else {
                // neither axis is skipped - split off first definition factor and compute multr
                if(termOper.isStandard(iIndex,jIndex)){
                    const BasisMatMatrix* mat=BasisMatMatrix::factory(termOper.first(),iIndex->basis(),jIndex->basis());

                    if(mat){
                        mult=mat->mat();
                    }
                    else if(BasisMatDependent::modify(termOper,iIndex,jIndex)!=termOper){
                        //                            mult=UseMatrix::Constant(iIndex->childSize(),jIndex->childSize(),1.);
                        mult=Eigen::MatrixXcd::Constant(iIndex->childSize(),jIndex->childSize(),1.);
                        termOper=BasisMatDependent::modify(termOper,iIndex,jIndex); // replace with modified operator
                        if(iIndex->axisName()=="Eta" and termOper.find("CO2")!=string::npos){
                            const BasisExpIm * b=dynamic_cast<const BasisExpIm*>(iIndex->parent()->basis());
                            if(b==0)ABORT("Hartree operator only for BasisExpIm (expIm for Phi-axis)");
                            int lim=20;
                            PrintOutput::DEVwarning(Sstr+"constrained to l<"+lim+"in"+termOper,1);
                            lim-=b->mValueAbs(iIndex->nSibling());
                            for(int j=0;j<mult.cols();j++)
                                for(int i=0;i<mult.rows();i++)
                                    if(i>lim or j>lim)mult(i,j)=0.;
                        }
                    }
                    else {
                        DEVABORT("cannot determine <"+idx()->strNode()+" | "+termOper+" | "+jdx()->strNode()+">");
                    }
                }
                else {
                    // we cannot make any guesses about the locality
                    int iSiz=iIndex->childSize(),jSiz=jIndex->childSize();
                    if(iIndex->hasFloor())iSiz=1;
                    if(jIndex->hasFloor())jSiz=1;
                    mult=Eigen::MatrixXcd::Constant(iSiz,jSiz,1.);
                }

                // descend in index on either side (if not floor)
                if(not iSub0->hasFloor())iSub0=iSub0->child(0);
                if(not jSub0->hasFloor())jSub0=jSub0->child(0);

                // constrain mult as given in definition (speeds up the setup)
                UseMatrix umult(mult);
                opSub=termOper.constrain(umult,iIndex,jIndex);
                mult=Eigen::Map<Eigen::MatrixXcd>(umult.data(),umult.rows(),umult.cols());


            }
            if(mult.size()==0)DEVABORT(Str("mult was left undefined in")+name+iIndex->strNode()+jIndex->strNode());
            if(iSub0->parent()->childSize()<mult.rows()+iSub0->nSibling())DEVABORT("mult does not match i");
            if(jSub0->parent()->childSize()<mult.cols()+jSub0->nSibling()){
                DEVABORT("mult does not match j");
            }

            const Index* jSub=jSub0;
            for(int j=0;j<mult.cols();j++,jSub=jSub->rightSibling()){
                const Index* iSub=iSub0;
                for(int i=0;i<mult.rows();i++,iSub=iSub->rightSibling()){
                    //                        if(abs(mult(i,j).complex())<epsPurge)continue;
                    if(abs(mult(i,j))<epsPurge)continue;
                    //                        childAdd(new OperatorTree(Name,opSub,iSub,jSub,Multiplier*mult(i,j),TFac,this));
                    addMulti(Name,opSub,iSub,jSub,Multiplier*mult(i,j),TFac);
                }
            }

        }

    }
    Timer::generalMonitor.monitor(0,Sstr+name+" at "+iIndex->index());
    STOPDEBUG(opAll);
    Index::build=false;
}

static std::string plotAction(const OperatorAbstract & Op, const Coefficients * C=0, std::shared_ptr<Plot> Plt=0){
    if(Plt==0)Plt.reset(new Plot(Op.idx(),ReadInput::main));
    if(Plt->isEmpty())return "cannot plot action of "+Op.name+", neet Plot: ... in input";

    Coefficients c(Op.jdx()),d(Op.idx()),p(Op.idx());
    std::string mess=Op.name;
    const OperatorAbstract* op=&Op;
    Coefficients *b(0);
    if(C!=0)
        c=*C;
    else {
        // can only set to one on numerical discretiztion
        c.setToZero();
        b=&c;
        while(b and (b->idx()->isHybrid() or dynamic_cast<const BasisOrbital*>(b->idx()->basis())!=0))
            b=b->nodeNext();
        if(b==0)return "cannot plot "+Op.name+"on hierarchy"+C->idx()->hierarchy();
        b->setToFunction("1");
    }

    if(b==&c)
        Op.apply(1.,c,0.,d);
    else {
        const OperatorTree* ot=dynamic_cast<const OperatorTree*>(op);
        mess=Op.name+"on "+b->idx()->hierarchy()+" of hierarchy"+C->idx()->hierarchy();
        while(ot and ot->jdx()!=b->idx() and op->idx()!=b->idx())ot=ot->nodeNext();
        if(ot==0)return "cannot apply"+mess;
        ot->apply(1,*b,0,*d.retrieve(ot->idx()));
    }
    Op.idx()->inverseOverlap()->apply(1,d,0.,p);

    std::string file="plot_"+Op.name+(C==0?"":"_C");
    file=ReadInput::main.output()+file;
    Plt->withPlotReal().plot(p,file,{"# "+Op.name+" "+Op.def().str()});
    return "action of "+mess+" on "+file;
}

void OperatorTree::postProcess(){

    // find (all) channel nodes
    for (std::string kind: {"XC","Hartree"}){
        if(definition.str().find(kind)!=std::string::npos){
            LOG_PUSH("postProcess");
            for(OperatorTree* op=this;op!=0;op=op->nodeNext()){
                //                if(OperatorData::terms(op->definition).size()==1
                if(OperatorDefinition(op->definition).terms().size()==1
                        and op->definition.str().find(kind)!=string::npos
                        and (op->iIndex->axisName()=="Channel" or op->iIndex->axisName()=="Orbital&Phi")
                        and op->iIndex->basis()==op->jIndex->basis()
                        ){

                    std::string pot=definition.str().substr(definition.str().find(kind));
                    pot=pot.substr(0,pot.find(">"));
                    pot=tools::stringInBetween(pot,"[","]");
                    if(pot==kind)pot="CoulombEE"; // default potential

                    std::vector<std::vector<Eigen::MatrixXcd>> rho;
                    const BasisOrbitalNumerical* orbs(0);
                    if(dynamic_cast<const BasisChannel*>(op->iIndex->basis())){
                        const BasisChannel* bI=dynamic_cast<const BasisChannel*>(op->iIndex->basis());
                        rho=bI->rho1();
                        orbs=dynamic_cast<const BasisOrbitalNumerical*>(bI->orbs());
                    }

                    if(kind=="XC")OperatorFloorXC::postProcess(this,pot,rho,orbs);
                    else if(kind=="Hartree")OperatorHartree::postProcess(this,rho,orbs);
                    else DEVABORT("define post-processing for kind="+kind);
                }
            }
            LOG_POP();
        }
    }

    LOG("floor");

    STARTDEBUG(operPar);
    ParallelOperator par(this);
    STOPDEBUG(operPar);
    STARTDEBUG(fuseOptree);
    par.fuse();
    STOPDEBUG(fuseOptree);


    STARTDEBUG(purgeOptree);
    par.purge();
    STOPDEBUG(purgeOptree);

    LOG("floor");

    //    if(name.find("Exp")==0 and (def().find("Hartree")!=string::npos or def().find("1*Q*Q")!=string::npos)){
    if(name.find("Exp")==0 and def().find("Hartree")!=std::string::npos){

        const Index* ix=idx();
        if(ix->axisName().find("&")!=std::string::npos)ix=ix->child(1);

        const OperatorTree* op=this;
        while(op and (op->idx()!=ix or op->jdx()!=ix))op=op->nodeNext();
        if(op)PrintOutput::message(plotAction(*op));
        else PrintOutput::DEVwarning("could not plot "+name+" potentials, no plotable block found");

    }


    MPIwrapper::Barrier();
    par.setDistribution();
    PrintOutput::progressStop();
}

static string fail;
void OperatorTree::fuse(){

    if(not OperatorAbstract::fuseOp)return; // globally switched off

    // vector<double>maxNorm;
    for (unsigned int k=0;k<childSize();k++){
        if(child(k)==0)ABORT("this should not happen");
        // maxNorm.push_back(child(k)->norm());
        for(unsigned int l=childSize()-1;l>k;l--){
            // maxNorm[k]=max(maxNorm[k],child(l)->norm()); // keep track of norms before fusing
            if(child(k)->absorb(child(l)))childEraseNode(l);
        }
    }
}
void OperatorTree::fuseBottomUp(){
    if(not OperatorAbstract::fuseOp)return; // globally switch off
    for(size_t k=0;k<childSize();k++)child(k)->fuseBottomUp();
    fuse();
}

bool OperatorTree::absorb(OperatorTree *Other){
    // do not absorb on ion level: possibly basis Ion has incorrect index.

    // only operator with equal indices, equal time-dependence, and equal tensor can be absorbed
    if(iIndex!=Other->iIndex)return false;
    if(jIndex!=Other->jIndex)return false;
    if(oFloor!=0){
        if((iIndex->isHybrid() || jIndex->isHybrid()) && (oFloor->kind()=="DUM")!=(Other->oFloor->kind()=="DUM"))
            PrintOutput::DEVwarning(Str("absorb")+"\n"+strNode(-1)+iIndex->basis()->str()
                                    +"\n"+Other->strNode(-1)+Other->iIndex->basis()->str()+"\n");
        if(Other->iIndex->root()->hierarchy().find("ValDer") != string::npos)return false;
        if(not OperatorFloor::absorb(oFloor,Other->oFloor,Str(definition,"")+Other->definition+iIndex->label()+jIndex->label()))return false;
    }

    if(name.find("(fused)")==string::npos)name+="(fused)";
    definition+=" "+Other->definition;

    for(unsigned int l=0;l<Other->childSize();l++)childAdd(Other->child(l));

    fuse();
    return true;
}

void OperatorTree::floorInvert(){

    for(OperatorTree* oLeaf=firstLeaf();oLeaf!=0;oLeaf=oLeaf->nextLeaf()){
        if (oLeaf->oFloor!=0 && dynamic_cast<OperatorDUM*>(oLeaf->oFloor) == 0){
            if(oLeaf->oFloor!=0){
                UseMatrix mat;
                oLeaf->oFloor->matrix(mat);
                if(mat.rows()!=mat.cols())ABORT(Sstr+"cannot invert floor"+oLeaf->index()+"of"+name+"at"+iIndex->strNode()+jIndex->strNode());
                mat=mat.inverse();
                delete oLeaf->oFloor;
                oLeaf->oFloor=OperatorFloor::factory(
                            std::vector<const UseMatrix*>(1,const_cast<const UseMatrix*>(&mat)),
                            Str(oLeaf->def())+"_inv:"+this);
            }
        }
    }
}

OperatorTree & OperatorTree::addColumnwise(OperatorTree *Block, bool View){

    std::vector<const Index*> jpath=Block->jIndex->path();
    std::vector<const Index*> ipath=Block->iIndex->path();

    if(Block->iIndex->root()!=iIndex)DEVABORT("lhs index roots do not match");
    if(Block->jIndex->root()!=jIndex)DEVABORT("rhs index roots do not match");

    OperatorTree* node=this;

    for(size_t p=0;p<jpath.size();p++){
        size_t k=0;
        while(k<node->childSize() and node->child(k)->jIndex->nSibling()<Block->jIndex->nSibling())k++;
        if(k==node->childSize() or node->child(k)->jIndex!=Block->jIndex){
            node->childInsert(k,new OperatorTree("n/a",ipath[0],jpath[p]));
        }
        node=node->child(k);
    }
    for(size_t p=0;p<ipath.size();p++){
        size_t k=0;
        while(k<node->childSize() and node->child(k)->iIndex->nSibling()<Block->iIndex->nSibling())k++;
        if(k==node->childSize() or node->child(k)->iIndex!=Block->iIndex)
            node->childInsert(k,new OperatorTree("n/a",ipath[p],jpath.back()));
        node=node->child(k);
    }
    if(View)node->childView(Block);
    else    node->childAdd(Block);
    return *node->childBack();
}


OperatorTree & OperatorTree::add(const OperatorTree *Term){

    // nothing to be added
    if(Term->isZero())
        goto Return;

    // add to zero
    if(isZero()){
        for(size_t k=0;k<Term->childSize();k++){
            childAdd(Term->child(k)->deepCopy());
            goto Return;
        }
    }

    if(iIndex!=Term->iIndex or jIndex!=Term->jIndex)
        ABORT("adding requires identical index blocks");

    if(Term->isLeaf() or isLeaf()){
        cout<<"Term\n"<<Term->str()<<endl;
        cout<<"leftn"<<str()<<endl;
        ABORT("cannot add leaf's");
    }

    if(!iIndex->isHybrid()) {
        if(Term->height()>height()){
            if(Term->child(0)->iIndex!=Term->iIndex or Term->child(0)->jIndex!=Term->jIndex){
                ABORT("hierarchies differ");
            }
            // insert dummy hiearchy level in present
            childAdd(new OperatorTree(name,iIndex,jIndex));
            childBack()->nodeCopy(this,false);
            for(size_t k=0;k<childSize()-1;k++)
                childBack()->childAdd(child(k));
            for(int k=childSize()-2;k>=0;k--)childEraseNode(k);
        }
        if(Term->height()!=height()){
            cout<<"Term1"<<str()<<endl;
            cout<<"Term2"<<Term->str()<<endl;
            ABORT("adding requires identical operator hierarchies (for now)");
        }
    }

    name+="+"+Term->name;
    definition+=" + ("+Term->definition+")";
    for(size_t k=0;k<Term->childSize();k++){
        childAdd(Term->child(k)->deepCopy());
        childBack()->name += Term->child(k)->name;
        childBack()->definition += Term->child(k)->definition;
    }
Return:
    ParallelOperator par(this);
    par.fuse();
    par.purge();
    return *this;
}

void OperatorTree::nodeCopy(const OperatorTree* Other, bool View){
    iIndex = Other->iIndex;
    jIndex = Other->jIndex;
    oFloor = Other->oFloor;
    if(not View and Other->oFloor!=0){
        if     (dynamic_cast<OperatorZG*>(Other->oFloor)!=0)oFloor=new OperatorZG(*dynamic_cast<OperatorZG*>(Other->oFloor));
        else if(dynamic_cast<OperatorZD*>(Other->oFloor)!=0)oFloor=new OperatorZD(*dynamic_cast<OperatorZD*>(Other->oFloor));
        else if(dynamic_cast<OperatorZGxZG*>(Other->oFloor)!=0)oFloor=new OperatorZGxZG(*dynamic_cast<OperatorZGxZG*>(Other->oFloor));
        else if(dynamic_cast<OperatorZDxZG*>(Other->oFloor)!=0)oFloor=new OperatorZDxZG(*dynamic_cast<OperatorZDxZG*>(Other->oFloor));
        else if(dynamic_cast<OperatorZGxZD*>(Other->oFloor)!=0)oFloor=new OperatorZGxZD(*dynamic_cast<OperatorZGxZD*>(Other->oFloor));
        else if(dynamic_cast<OperatorZero*>(Other->oFloor)!=0)oFloor=new OperatorZero(*dynamic_cast<OperatorZero*>(Other->oFloor));
        else if(dynamic_cast<OperatorFloorEE*>(Other->oFloor)!=0)oFloor=new OperatorFloorEE(*dynamic_cast<OperatorFloorEE*>(Other->oFloor));
        else if(dynamic_cast<OperatorZGSblock*>(Other->oFloor)!=0)oFloor=new OperatorZGSblock(*dynamic_cast<OperatorZGSblock*>(Other->oFloor));
        else if(dynamic_cast<OperatorDUM*>(Other->oFloor)!=0)oFloor=new OperatorDUM(*dynamic_cast<OperatorDUM*>(Other->oFloor));
        else
            PrintOutput::DEVwarning("in OperatorTree::nodeCopy -- "
                                    +Other->oFloor->kind()
                                    +" not copied, only pointed to -- may break, if one of the terms is deleted",1);
    }
}

#include <typeinfo>
bool OperatorTree::nodeEquivalent(const OperatorTree* Other) const{
    if(iIndex!=Other->iIndex)return false;
    if(jIndex!=Other->jIndex)return false;
    if(oFloor==Other->oFloor)return true;
    if((!oFloor and Other->oFloor) or (oFloor and !Other->oFloor))return false;
    if(oFloor->kind()!=Other->oFloor->kind())return false;
    if(oFloor->matrix()==Other->oFloor->matrix())return true;
    return false;
}

//OperatorTree::OperatorTree(const Operator* A)
//    :OperatorAbstract(A->name,A->iIndex,A->jIndex),oFloor(A->oFloor),_view(true)
//{
//    definition=A->definition;
//    // handle legacy operator
//    if(A->o.size()>0){
//        _view=false;
//        UseMatrix mat;
//        A->o[0]->matrix(mat);
//        oFloor=OperatorFloor::factory(vector<const UseMatrix*>(1,&mat),Str("legacyOverlap")+iIndex->label()+"_"+jIndex->label());
//    }
//    if(A->tensor!=0 or A->preTensor!=0)ABORT("tensor structured operator cannot be converted to OperatorTree (yet)");
//    for(size_t k=0;k<A->O.size();k++)childAdd(new OperatorTree(A->O[k]));
//}

UseMatrix OperatorTree::matrixBlocks(unsigned int IDepth, unsigned int JDepth) const {
    // determine matrix sizes
    int idim=1,jdim=1;
    // input depth is relative to presen node, internal depth from top of tree
    IDepth+=iIndex->depth();
    JDepth+=jIndex->depth();
    const Index *idx=iIndex,*jdx=jIndex;
    while(idx->descend()!=0 and not idx->hasFloor() and idx->depth()<IDepth)idx=idx->descend();
    while(jdx->descend()!=0 and not jdx->hasFloor() and jdx->depth()<JDepth)jdx=jdx->descend();
    IDepth=idx->depth();
    JDepth=jdx->depth();
    while(0!=(idx=idx->nodeRight(iIndex)))idim++;
    while(0!=(jdx=jdx->nodeRight(jIndex)))jdim++;

    const OperatorTree* opp=this;
    while(opp->iIndex->depth()<IDepth or opp->jIndex->depth()<JDepth)opp=opp->descend();

    UseMatrix mat=UseMatrix::Constant(idim,jdim,0.);
    for(const OperatorTree* b=this;b!=0;b=b->nodeNext(this))
        if(b->iIndex->depth()==IDepth and b->jIndex->depth()==JDepth){
            int i=b->iIndex->levelRank(iIndex->parent()),j=b->jIndex->levelRank(jIndex->parent());
            mat(i,j)=max(mat(i,j).real(),b->norm());
        }
    return mat;
}

void OperatorTree::matrix(UseMatrix &Mat) const{
    if(isHuge())ABORT("too big - cannot construct full matrix");
    if(Mat.size()==0)
        Mat=UseMatrix::Zero(iIndex->sizeStored(),jIndex->sizeStored());
    if(isLeaf())
        OperatorAbstract::matrix(Mat);
    else{
        for(size_t k=0;k<childSize();k++){
            UseMatrix mat;
            child(k)->matrix(mat);
            Mat.block(child(k)->iIndex->posIndex(iIndex),child(k)->jIndex->posIndex(jIndex),mat.rows(),mat.cols())+=mat;
        }
    }
}
Eigen::MatrixXcd OperatorTree::matrix() const {
    Eigen::MatrixXcd mat;
    return matrix(mat,iIndex,jIndex);
}

Eigen::SparseMatrix<std::complex<double>> OperatorTree::matrixSparse(bool Contract, double Eps) const {
    if(childSize()==0 and isHuge())ABORT("huge (sub-)matrix");
    BlockView bv(this);
    Eigen::SparseMatrix<std::complex<double>> res;
    bv.sparseMatrix(res,Contract,0.);
    return res;
}

Eigen::MatrixXcd  & OperatorTree::matrix(Eigen::MatrixXcd & Mat, const Index* ISub, const Index* JSub) const {
    if(isHuge())ABORT("too big - cannot construct full matrix");
    if(ISub==0)ISub=iIndex;
    if(JSub==0)JSub=jIndex;
    if(Mat.rows()!=ISub->sizeStored() and Mat.cols()!=JSub->sizeStored()){
        if(Mat.size()!=0)DEVABORT("entering with non-matching matrix");
        Mat=Eigen::MatrixXcd::Zero(ISub->sizeStored(),JSub->sizeStored());
    }
    if(oFloor!=0){
        Eigen::MatrixXcd mat=oFloor->matrix();
        Mat.block(iIndex->posIndex(ISub),jIndex->posIndex(JSub),mat.rows(),mat.cols())+=mat;
    }
    else if (not childSize()){
        Mat.block(iIndex->posIndex(ISub),jIndex->posIndex(JSub),iIndex->size(),jIndex->size())+=
                OperatorAbstract::matrix();
    }
    else {
        for(size_t k=0;k<childSize();k++){
            // descend if sub-index is in current or current is part of subtree
            if(     (child(k)->iIndex->isSubtree(ISub) or ISub->isSubtree(child(k)->iIndex) ) and
                    (child(k)->jIndex->isSubtree(JSub) or JSub->isSubtree(child(k)->jIndex) ) )
                child(k)->matrix(Mat,ISub,JSub);
            else {
            }
        }
    }
    return Mat;
}

OperatorTree::OperatorTree(std::string Name, const Index* Idx, const Index* Jdx, const Eigen::MatrixXcd& Mat)
    :OperatorTree(Name,Idx,Jdx)
{
    _buildFromMatrix(Idx,Jdx,Mat);
}

void OperatorTree::_buildFromMatrix(const Index *Idx, const Index *Jdx, const Eigen::MatrixXcd &Mat){
    if(Idx->size()!=Mat.rows() or Jdx->size()!=Mat.cols())
        DEVABORT(Sstr+"matrix dimension does not match indices"
                 +Mat.rows()+" x "+Mat.cols()+" != "+Idx->size()+" x "+Jdx->size());
    if(Idx->hasFloor() and Jdx->hasFloor()){
        // !! caution - factory may rely on contigous storage of data !!
        oFloor=OperatorFloor::factory({&Mat},name+Idx->hash()+"|"+Jdx->hash());
    }
    else if(not Idx->hasFloor())
        for(size_t k=0;k<Idx->childSize();k++){
            childAdd(new OperatorTree(name,Idx->child(k),Jdx,
                                      Mat.block(Idx->child(k)->posIndex(Idx),0,Idx->child(k)->size(),Mat.cols())));
        }
    else {
        for(size_t k=0;k<Jdx->childSize();k++){
            childAdd(new OperatorTree(name,Idx,Jdx->child(k),
                                      Mat.block(0,Mat.rows(),Jdx->child(k)->posIndex(Jdx),Jdx->child(k)->size())));
        }
    }
}

void OperatorTree::reFloor(size_t IFloor, size_t JFloor){
    if(idx()->depth()==IFloor and jdx()->depth()==JFloor and floor()){
        if(childSize())DEVABORT("operator has floor and children?");
        if(not (idx()->hasFloor() and jdx()->hasFloor()))DEVABORT("operator floor, but not index floor?");
        return;
    }

    if(idx()->depth()<IFloor and idx()->depth()<JFloor)
        for(size_t k=0;k<childSize();k++)child(k)->reFloor(IFloor,JFloor);

    Eigen::MatrixXcd mat=matrix();
    const_cast<Index*>(idx())->resetFloor(IFloor);
    const_cast<Index*>(jdx())->resetFloor(JFloor);

    // erase all children
    for(size_t k=childSize();k>0;k--)childErase(k-1);

    // rebuild operator from matrix
    _buildFromMatrix(idx(),jdx(),mat);
}


void OperatorTree::_blockContracted(UseMatrix &Mat, int &I0, int &J0, const Index* IRoot, const Index* JRoot,
                                    const std::vector<unsigned int> &ICont, const std::vector<unsigned int> &JCont) const{
    int ipos=iIndex->posIndex(IRoot);
    int jpos=jIndex->posIndex(JRoot);
    // get block position in contracted matrix
    I0=*std::min_element(ICont.begin()+ipos,ICont.begin()+ipos+iIndex->size());
    J0=*std::min_element(JCont.begin()+jpos,JCont.begin()+jpos+jIndex->size());

    UseMatrix mat,cmat;
    matrix(mat);
    matrixContract(mat,cmat);
}

void OperatorTree::_matrixContractedOLD(UseMatrix &Mat) const{
    // entry level - define the numbering and create matrix
    vector<unsigned int> iCont=iIndex->contractedNumbering();
    vector<unsigned int> jCont=jIndex->contractedNumbering();
    // last index is largest
    Mat=UseMatrix::Zero(iCont.back()+1,jCont.back()+1);
    _matrixContracted(Mat,iIndex,jIndex,iCont,jCont);
}

Eigen::MatrixXcd OperatorTree::matrixContracted() const {
    UseMatrix mat;
    //    _matrixContractedOLD(mat);
    // entry level - define the numbering and create matrix
    std::vector<unsigned int> iCont=iIndex->contractedNumbering();
    std::vector<unsigned int> jCont=jIndex->contractedNumbering();
    // last index is largest - determines dimensions
    mat=UseMatrix::Zero(iCont.back()+1,jCont.back()+1);
    _matrixContracted(mat,iIndex,jIndex,iCont,jCont);
    return Eigen::Map<Eigen::MatrixXcd>(mat.data(),mat.rows(),mat.cols());
}


TIMER(contractM,)
TIMER(cMat,)
void OperatorTree::_matrixContracted(UseMatrix &Mat, const Index *IRoot, const Index *JRoot, const std::vector<unsigned int> &ICont, const std::vector<unsigned int> &JCont) const
{
    if(iIndex->continuity()!=Index::npos or jIndex->continuity()!=Index::npos){
        // continuity()!=Index::npos indicates that a level is to be contracted
        UseMatrix cmat;
        int i0,j0;
        _blockContracted(cmat,i0,j0,IRoot,JRoot,ICont,JCont);
        Mat.block(i0,j0,cmat.rows(),cmat.cols())=cmat;
    }
    else
        for(size_t k=0;k<childSize();k++)
            child(k)->_matrixContracted(Mat,IRoot,JRoot,ICont,JCont);
}

OperatorTree::OperatorTree(const OperatorAbstract *A, const Index* IIndex , const Index* JIndex)
    :OperatorAbstract(A->name,IIndex,JIndex),oFloor(0),_view(false)
{

    if(iIndex==0)iIndex=A->iIndex;
    if(jIndex==0)jIndex=A->jIndex;

    if(not jIndex->hasFloor()){
        for(size_t k=0;k<jIndex->childSize();k++){
            childAdd(new OperatorTree(A,iIndex,jIndex->child(k)));
            if(childBack()->isZero())childPop();
        }
    }

    else {
        if(iIndex==A->iIndex){
            // get matrix columns
            //            vector<complex<double> > colStor;
            colStor.clear();
            A->subMatrix(colStor,iIndex,jIndex);
            colMat=new UseMatrix(UseMatrix::UseMap(colStor.data(),iIndex->sizeStored(),jIndex->sizeStored()));
        }
        if(iIndex->hasFloor()){
            // place block into floor
            UseMatrix fMat;
            fMat=colMat->block(iIndex->posIndex(A->iIndex),0,iIndex->sizeStored(),jIndex->sizeStored());
            oFloor=OperatorFloor::factory(vector<const UseMatrix*>(1,&fMat.purge()),A->name);
        }
        else {
            for(size_t k=0;k<iIndex->childSize();k++){
                childAdd(new OperatorTree(A,iIndex->child(k),jIndex));
                if(childBack()->isZero())childPop();
            }
        }
        if(iIndex==A->iIndex){
            delete colMat;
        }
    }

}

void OperatorTree::write(std::ofstream &File){
    if(not File.is_open())DEVABORT("cannot write - file not open");

    // write name, definition, and indices
    if(File.tellp()==0){
        if(MPIwrapper::isMaster()){
            OperatorAbstract::write(&File);

            // write list of definitions of parameter-depedent factor
            tools::write(File,int(0));
            tools::write(File,functionStrings.size());
            for(auto fac: functionStrings){
                tools::write(File,fac.first);
                tools::write(File,fac.second.size());
                tools::write(File,fac.second.data(),fac.second.size());
            }
        }
    }

    if(floor()){
        int host=dynamic_cast<OperatorDUM*>(floor())?-1:MPIwrapper::Rank();
        MPIwrapper::AllreduceMAX(&host,1);
        // determine host
        std::vector<int> info(5,0);
        std::vector<std::complex<double>> buf;
        if(MPIwrapper::isMaster()){
            tools::write(File,int(2));
            OperatorFloor* f=floor();
            if(host!=MPIwrapper::master()){
                MPIwrapper::Recv(info.data(),info.size(),host,1);
                buf.resize(info[4]);
                MPIwrapper::Recv(buf.data(),buf.size(),host,2);
                f=OperatorFloor::unpackFactory(info,buf);
            }
            f->write(File);
            tools::write(File,floor()->factor());
            if(host!=MPIwrapper::master())delete f;
        }
        else if(MPIwrapper::Rank()==host) {
            floor()->pack(info,buf);
            MPIwrapper::Send(info.data(),info.size(),MPIwrapper::master(),1);
            MPIwrapper::Send(buf.data(),buf.size(),MPIwrapper::master(),2);
        }
        MPIwrapper::Barrier();

    }
    else {
        if(MPIwrapper::isMaster()){
            tools::write(File,int(1));
            tools::write(File,int(childSize()));
        }
        for(size_t k=0;k<childSize();k++){
            if(MPIwrapper::isMaster()){
                tools::write(File,int(child(k)->iIndex==iIndex?-1:child(k)->iIndex->nSibling()));
                tools::write(File,int(child(k)->jIndex==jIndex?-1:child(k)->jIndex->nSibling()));
            }
            child(k)->write(File);
        }
    }

}

static std::map<std::complex<double>*,std::complex<double>*> functionPtrs;
OperatorTree::OperatorTree(std::ifstream &File, string Name, const Index *IIndex, const Index *JIndex)
    :OperatorAbstract(File,Name,IIndex,JIndex),oFloor(0),_view(false)
{
    int code; // 0 - header, 1-operator, 2-floor
    tools::read(File,code);
    if (code==0){
        // get list of time dependent factors
        size_t n;
        tools::read(File,n);
        std::complex<double>* pFac;
        for(size_t k=0;k<n;k++){
            tools::read(File,pFac);
            size_t len;
            tools::read(File,len);
            std::string fString(len,' ');
            for(auto &c: fString)tools::read(File,&c,sizeof(c));

            functionPtrs[pFac]=Parameters::pointer(fString);
        }
        tools::read(File,code);
    }

    if(code==2){
        OperatorFloor* f=OperatorFloor::readFactory(File);
        std::complex<double>* pFac;
        tools::read(File,pFac);

        oFloor=Parallel::operatorFloor(IIndex,JIndex,[&](){ return f;});
        oFloor->setFactor(functionPtrs[pFac]);
        if(dynamic_cast<OperatorDUM*>(oFloor))delete f;
    }
    else if (code==1){
        int siz,ic,jc;
        tools::read(File,siz);
        for(int k=0;k<siz;k++){
            tools::read(File,ic);
            tools::read(File,jc);
            childAdd(new OperatorTree(File,name,
                                      ic==-1?iIndex:iIndex->child(ic),
                                      jc==-1?jIndex:jIndex->child(jc)));
        }
    }
    else {
        DEVABORT("illegal code");
    }
}

void OperatorTree::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const {
    if(&Vec==&Y)DEVABORT("Vec and Y must not be the same vector");
    int host=ParallelOperator::getHost(this);

    if(B!=0.){
        bool crit=timeCritical::coefficients;
        timeCritical::coefficients=false;
        *tempLHS()=Y;
        timeCritical::coefficients=crit;
    }
    // apply only if overall, distributed or local to thread
    if(        host==ParallelOperator::distributed
               or host==ParallelOperator::all
               or host==MPIwrapper::Rank()
               ){
        _apply(A,Vec,0.,Y);
    }

    if(host==ParallelOperator::distributed and MPIwrapper::Size()>1){
        Y.treeOrderStorage();
        // HACK: in Hybrid Disc, if Y is Ionic coefficient, it will still yield Y.data() == 0, because
        // the contiguous floor data is put at the Hybrid hierarchy level
        // --> needs more thoughtful fix!
        if(Y.data() == 0){
            MPIwrapper::AllreduceSUM(const_cast<Coefficients*>(Y.parent())->data() + 1,Y.parent()->size() - 1);
        } else {
            MPIwrapper::AllreduceSUM(Y.data(),Y.size());
        }
    }
    if(B!=0.)Y.axpy(B,tempLHS());
}

void OperatorTree::_apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    Y*=B;
    if(idx()!=Y.idx()){
        ABORT(Sstr+"lhs index does not match: oper"+idx()+"\n"+idx()->str()+"\noutput"+Y.idx()+"\n"+Y.idx()->str(0));
    }
    if(A==0. or Vec.isZero())return;

    if(jIndex!=Vec.idx() and not jIndex->hasFloor())DEVABORT("rhs index does not match: \nOper\n"+jIndex->str(0)+"\nVec\n"+Vec.idx()->str(0));
    if(floor()!=0 and floor()->kind()!="DUM"){
        oFloor->apply(A,Vec.data(),Vec.size(),1.,Y.data(),Y.size());
        return;
    }
    else if(not childSize())
        // a leaf w/o a floor may by a derived class operator
        applyDerivedClass(A,Vec,B,Y);
    else
        for (size_t n=0;n<childSize();n++){
            const Coefficients* jVec=&Vec;
            Coefficients* iVec = &Y;

            if(child(n)->iIndex!=iIndex)iVec=iVec->child(child(n)->iIndex->nSibling());
            if(child(n)->jIndex!=jIndex and not jIndex->hasFloor())jVec=jVec->child(child(n)->jIndex->nSibling());

            child(n)->_apply(A,*jVec,1.,*iVec);
        }
}

void OperatorTree::applyDerivedClass(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
{
    if(not isLeaf())DEVABORT("cannot applyDerivedClass for non-leaf");
    if(dynamic_cast<const OperatorDiagsPermuted*>(this))
        dynamic_cast<const OperatorDiagsPermuted*>(this)->apply(A,Vec,B,Y);
    if(dynamic_cast<const OperatorPermute*>(this))
        dynamic_cast<const OperatorPermute*>(this)->apply(A,Vec,B,Y);
}

void OperatorTree::applyTranspose(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    Y*=B;
    if(jIndex!=Y.idx())ABORT("rhs index does not match input Vector");
    if(iIndex!=Vec.idx())ABORT("lhs index does not match input Vector");
    if(oFloor!=0){
        Coefficients* jVec=const_cast<Coefficients*>(&Vec); // not nice...
        vector<complex<double> > vec(jVec->data(),jVec->data()+Vec.size());
        vector<complex<double> > y(Y.data(),Y.data()+Y.size());
        ABORT("re-implement this");
        //        oFloor->axpyTranspose(A,vec,1.,y);
        for(size_t k=0;k<Y.size();k++)Y.data()[k]=y[k];
        return;
    }

    for (size_t n=0;n<childSize();n++){
        Coefficients* jVec=&Y;
        const Coefficients* iVec = &Vec;
        if(child(n)->iIndex!=iIndex)iVec=iVec->child(child(n)->iIndex->nSibling());
        if(child(n)->jIndex!=jIndex)jVec=jVec->child(child(n)->jIndex->nSibling());
        child(n)->applyTranspose(A,*iVec,1.,*jVec);
    }
}

bool OperatorTree::isSymmetric(std::string Kind, double Eps) const {
    BlockView b(this);
    return b.isSymmetric(Kind,Eps);
}

bool OperatorTree::isDiagonal() const {

    if(isLeaf())
        return floor()==0 or (iIndex==jIndex and floor()->isDiagonal());

    for(size_t k=0;k<childSize();k++)
        if(child(k)->iIndex!=child(k)->jIndex or not child(k)->isDiagonal())
            return false;

    return true;
}

bool OperatorTree::isBlockDiagonal() const {
    if(isLeaf())return false;
    for(size_t k=0;k<childSize();k++){
        if(child(k)->idx()!=child(k)->jdx())return false;
        if(child(k)->childSize()>1 and child(k)->descend()->idx()==child(k)->idx())return false; // multi-component cannot be handled
    }
    return iIndex->continuity()==Index::npos;
}

bool OperatorTree::isZero(double Eps, bool Stochastic) const {
    if(iIndex->sizeStored()==0)return true;
    if(jIndex->sizeStored()==0)return true;
    if(not childSize()){
        if(not oFloor){
            bool res=OperatorAbstract::isZero(Eps,Stochastic);
            return res;
        }
        else
            return oFloor->norm()<=Eps;
    }
    for(size_t k=0;k<childSize();k++)
        if(not child(k)->isZero(Eps,Stochastic))return false;
    return false;
}

Coefficients  OperatorTree::diagonal(bool OnlyForDiagonalOperator) const{
    if(OnlyForDiagonalOperator and iIndex!=jIndex)ABORT("operator is not  block diagonal");

    Coefficients res(iIndex,0.);
    if(childSize()>0 and child(0)->iIndex==iIndex){
        // in case there are multiple diagonal blocks...
        res=child(0)->diagonal(OnlyForDiagonalOperator);
        for(size_t k=1;k<childSize();k++)res+=child(0)->diagonal(OnlyForDiagonalOperator);
    }

    else {
        if(oFloor!=0){
            if(floor()->cols()>0){
                // place diagonal into coefficients
                UseMatrix mat;
                oFloor->matrix(mat);

                if(OnlyForDiagonalOperator and not mat.isDiagonal(1.e-12)){
                    mat.print(Str("mat\n")+iIndex->str()+jIndex->str(),2);
                    ABORT("operator floor is not diagonal");
                }
                // not fast, but save
                for(size_t k=0;k<min(iIndex->sizeStored(),jIndex->sizeStored());k++)
                    res.data()[k]=mat(k,k).complex();
            }
            if(ParallelOperator::getHost(this)==ParallelOperator::distributed)
                MPIwrapper::AllreduceSUM(res.data(),res.size());
        }
        for(size_t k=0;k<childSize();k++){
            if(child(k)->iIndex==child(k)->jIndex)
                *res.child(child(k)->iIndex->nSibling())+=child(k)->diagonal(OnlyForDiagonalOperator);
            else
                if(OnlyForDiagonalOperator)ABORT("operator is not block-diagonal");
        }
    }
    return res;
}

double OperatorTree::norm() const {
    if(iIndex->sizeStored()==0)return 0.;
    if(jIndex->sizeStored()==0)return 0.;

    double nrm;
    if(oFloor!=0)
        nrm=oFloor->norm();
    else {
        nrm=0.;
        for(size_t k=0;k<childSize();k++)nrm=max(nrm,child(k)->norm());
    }
    return nrm;
}

void OperatorTree::replaceIndex(const Index* IRep, const Index* JRep){

    if((not IRep or idx()==IRep) and (not JRep or jdx()==JRep))return;
    if(IRep and idx()->childSize()!=IRep->childSize())DEVABORT("lhs replacement does not match structure");
    if(JRep and jdx()->childSize()!=JRep->childSize())DEVABORT("rhs replacement does not match structure");


    for(size_t n=0;n< childSize();n++){
        const Index *iRep(IRep),*jRep(JRep);
        if(iRep and child(n)->idx()!= idx())iRep=iRep->child(child(n)->idx()->nSibling());
        if(jRep and child(n)->jdx()!= jdx())jRep=jRep->child(child(n)->jdx()->nSibling());
        child(n)->replaceIndex(iRep,jRep);
    }
    if(IRep)iIndex=IRep;
    if(JRep)jIndex=JRep;
}

void OperatorTree::purge(double Eps){
    for(int k=childSize();k>0;k--){
        child(k-1)->purge(Eps);
        if(child(k-1)->OperatorTree::isZero(Eps))childErase(k-1);
    }
}

void OperatorTree::purge(purgeCriterion Crit){
    for(int k=childSize();k>0;k--){
        child(k-1)->purge(Crit);
        if((*Crit)(child(k-1)))childErase(k-1);
    }
}

const OperatorFloor* OperatorTree::floor() const{
    return const_cast<const OperatorFloor*>(oFloor);
}

OperatorFloor*& OperatorTree::floor() {
    return oFloor;
}

std::string OperatorTree::strNode(int Digits) const {
    bool ptrs=Digits==Tree_withPtrs;
    if(ptrs)Digits=Tree_defaultKind;
    Str s("","");
    if(parent()==0)s=s+name+"\n";
    if(MPIwrapper::Size()>1)s=s+" <"+ParallelOperator::strHost(this)+">";
    if(isView())s=s+"V";
    if(ptrs)s+=Str("","")+iIndex+"<-"+jIndex+" ";
    s=s+tools::str(norm());
    s=s+" <"+iIndex->index()+"|"+jIndex->index()+"> ("+iIndex->sizeCompute()+"x"+jIndex->sizeCompute()+") ["+childSize()+"] "+iIndex->axisName();
    if(iIndex->axisName()!=jIndex->axisName())s+="-"+jIndex->axisName();
    //    if(iIndex->axisName()!=jIndex->axisName())s+="..."+iIndex->strNode()+"|"+jIndex->strNode();
    if(definition.length()>60)s+": "+definition.substr(0,60)+"...";
    else s=s+": "+definition;
    //    s=s+" ptr: "+iIndex+" | "+jIndex;
    if(oFloor!=0){
        Str proc("","");
        if(Digits==Tree_defaultKind or ptrs)s=s+" F("+oFloor->kind()+")"+proc;
        else                                s=s+proc+" (F) "+oFloor->str(Digits);
    }
    return std::move(s);
}

long OperatorTree::applyCount() const{
    if(oFloor!=0) return oFloor->applyCount();

    long result=0;
    for(unsigned int i=0; i<childSize(); i++){
        result+=child(i)->applyCount();
    }

    return result;
}

OperatorTree* OperatorTree::testOp(int Dim, string File){
    if(File=="")File=BasicDisc::generateTestInputs(0,0,Dim);
    const BasicDisc * D=new BasicDisc(File);

    string def="0.5<<Laplacian>>+0.5<<Parabolic>>";
    return new OperatorTree("testOp",def,D->idx(),D->idx());

}

size_t OperatorTree::diagnoseSizeOfNode() const {
    size_t sizAll=0;
    if(iIndex->parent()==0)
        PrintOutput::DEVmessage(Str("OperatorTree sizes: ")+memoryUsage["def"]+":"+sizAll);
    return sizAll;
}

void OperatorTree::update(double Time, const Coefficients *C){
    Parameters::update(Time);
}
void OperatorTree::updateNonLin(double Time, const Coefficients *C){
    if(floor()!=0 and floor()->kind()!="DUM"){
        OperatorFloorNonLin* fNL=dynamic_cast<OperatorFloorNonLin*>(floor());
        if(fNL)fNL->updateNonLin(_time, C);
        return;
    }
    for (size_t n=0;n<childSize();n++){
        child(n)->updateNonLin(Time, C);
    }
}


