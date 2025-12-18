// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include <set>

#include "tRecXchecks.h"
#include "discretizationHybrid.h"

#include "mpiWrapper.h"
#include "printOutput.h"
#include "readInput.h"
#include "basicDisc.h"
#include "index.h"
#include "operatorData.h"
#include "operatorVectors.h"
#include "inverse.h"
#include "operatorTree.h"
#include "basisNdim.h"
#include "basisDvr.h"

#ifdef _USE_HACC_
#include "discretizationHaCC.h"
#include "quantumChemicalInput.h"
#endif

#include "coefficientsGlobal.h"
#include "coefficientsLocal.h"
#include "parallel.h"
#include "parallelOperator.h"
#include "basisOrbital.h"
#include "axisTree.h"

using namespace std;

class BasicDiscDum: public BasicDisc{
    // dummy constructor for BasicDisc, avoid pollution of BasicDisc
public:
    BasicDiscDum(const Index* Idx){
        idx()=const_cast<Index*>(Idx);
    }
};


static std::set<std::vector<std::string> > knownHybrids=
{
    {"Neut","Chan"},
    {"Subspace","Complement"},
    {"Off","Center"},
    {"Added","Basic"}
};

DiscretizationHybrid::Sinverse::~Sinverse(){delete sInvBA; delete sAB; delete sAinv; delete sBinv; delete aVec; delete bVec;}

const Discretization* DiscretizationHybrid::haCC() const {
    if(comp[1]->getAxis()[0].name!="Ion")
        ABORT(Str("does not seem to be hybrid of haCC:")+idx()->hierarchy());
    return comp[1];
}

bool DiscretizationHybrid::isHybrid(ReadInput & Inp, int Line){
    AxisTree axT(Inp,Line);
    if(axT.height()<2)return false;
    return axT.child(0)->childSize()>1; // multiple axes on next level indicate hybrid
}



static std::unique_ptr<CoefficientsGlobal> gY, gX;

/** Example input:
 *>
 *> --- direct sum of discretizations named "Neut" and "Chan"\n
 *> Hybrid: subset=Neut Chan
 *>
 *> --- axes for the two components \n
 *> Axis: subset, name,nCoefficients,lower end, upper end,functions,order \n
 *> Neut, Neutral,1 \n
 *> Chan, X,20,-Infty,-20,polExp[1.] \n
 *>      ,X,10,-20, 20,polynomial,10 \n
 *>      ,X,20, 20, Infty,polExp[1.]
 *>
 *> --- operators must be specified block-wise for the components, here two diagonal blocks <0,0> and <1,1>:\n
 *> Operator: hamiltonian='-137<0,0><1>+<1,1>(<d_0.5_d>-<1/sqrt(Q*Q+2)>)'
 *>
 */
DiscretizationHybrid::DiscretizationHybrid(ReadInput &In)
{
    // construct components according to compName:
    if(not In.found("Axis","subset"))ABORT("specify 'Axis: subset' for using hybrid discretization");

    // collect hybrid components
    int line=0;
    string subs="";
    while(not In.endCategory("Axis",++line)){
        compName.push_back(AxisTree::readSubset(In,"Axis",line,subs));
        if(line>1 and compName.back()==subs)compName.pop_back();
        else subs=compName.back();
    }

    // get then names of the components
    std::string ovrOff="+1.<1,0>[[Id]]+1.<0,1>[[Id]]";
    if(find(knownHybrids.begin(),knownHybrids.end(),compName)==knownHybrids.end()){
        Str s("","\n");
        for(vector<string> c: knownHybrids)s+=c;
        ABORT(Str("Unknown hybrid components:","\n")+compName+"Available:"+s);
    }

    vector<string> cmp={"Subspace","Complement"};
    if(compName==cmp)ovrOff="";
    In.read("Hybrid","ovrOffdiag",ovrOff,ovrOff,"overlap blocks coupling the subsets");

#ifdef _USE_HACC_
    // special case treatment for haCC
    if(compName.size()==2 and compName[1]=="Neut" and compName[0]=="Chan")
        ABORT("Choose the other axis order: Neut Chan");
    if(compName[0]=="Neut")QuantumChemicalInput::read(In);
#endif

    // top level of index
    idx()=new Index();
    idx()->setAxisName("Hybrid");
    idx()->setAxisSubset(compName[0]+"&"+compName[1]);
    idx()->setBasis(BasisAbstract::factory("Hybrid: "+tools::str(compName.size())));

    // create the hybrid axis
    BasisSetDef bas;
    bas.coor=Coordinate::fromString("Hybrid");
    bas.order=compName.size();
    axis.push_back(Axis("",ComplexScaling(),bas));


    comp.resize(compName.size());
    for(int k=compName.size()-1;k>=0;k--){
        // replace below by factory
        if(compName[k]=="Chan"){
#ifdef _USE_HACC_
            DEVABORT("disabled for nowS");
//            comp[k]=new DiscretizationHaCC(In,compName[k],comp[k]->idx()->basisSet()->PointerToFunctions());
#else
            DEVABORT("for using DiscretizationHaCC, compile with -D_USE_HACC_");
#endif
        }
        else {
            comp[k]=new BasicDisc(In,compName[k]);
            BasisOrbital::referenceIndex[compName[k]]=comp[k]->idx();
        }
    }
    for(size_t k=0;k<compName.size();k++){
        idx()->childAdd(new Index(*comp[k]->idx()));
        idx()->childBack()->setAxisSubset(compName[k]);

        // compose names etc.
        name+= compName[k]+"."+comp[k]->name+" (+) ";
        axis[0].name+=compName[k]+"&";
        for(size_t m=0;m<comp[k]->getAxis().size();m++){
            axis.push_back(comp[k]->getAxis()[m]);
            axis.back().name+="("+compName[k]+")";
        }
    }

    name.resize(name.length()-5); // trim trailing " (+) "
    axis[0].name.resize(axis[0].name.length()-1);
    // get all index sizes
    idx()->sizeCompute();


    // top level of overlap operator
    // default overlap is direct sum: S = S_0 (+) S_1 (+)...(+) S_n-1
    idx()->localOverlapAndInverse(0,0);

    if(comp.size()>2)ABORT("for now only two components allowed");

    // add possible couplings between components
    //ovrOff="";
    Overlap * o=new Overlap(this,ovrOff);
    idx()->setOverlap(o);

    // get the inverse overlap oparator (block-inverse)
    if(not ReadInput::main.flag("DEBUGfem","run FEM basis - compute exact integrals everywhere"))
        idx()->setInverseOverlap(new Sinverse(this,o));
}

DiscretizationHybrid::DiscretizationHybrid(const Index* Idx)
{
    // construct components according to compName:
    //    if(not In.found("Axis","subset"))ABORT("specify 'Axis: subset' for using hybrid discretization");
    if(Idx->axisName().find("&")!=std::string::npos)DEVABORT("must have hybrid as first axis, got: "+Idx->axisName());

    idx()=const_cast<Index*>(Idx);
    for(size_t i=0;i<Idx->childSize();i++)compName.push_back(Idx->child(i)->axisSubset());
    if(compName.size()>2)ABORT("for now only two components allowed");

    // get then names of the components
    if(find(knownHybrids.begin(),knownHybrids.end(),compName)==knownHybrids.end()){
        Str s("","\n");
        for(vector<string> c: knownHybrids)s+=c;
        ABORT(Str("Unknown hybrid components:","\n")+compName+"Available:"+s);
    }

    vector<string> cmp={"Subspace","Complement"};
    if(compName!=cmp)DEVABORT(Sstr+"need "+cmp+", got:"+compName);
    //    In.read("Hybrid","ovrOffdiag",ovrOff,ovrOff,"overlap blocks coupling the subsets");

    for(int k=Idx->childSize()-1;k>=0;k--){
        BasisOrbital::referenceIndex[compName[k]]=Idx->child(k);
        comp.push_back(new BasicDiscDum(Idx->child(k)));
    }


    // add possible couplings between components
    std::string ovrOff="+1.<1,0>[[Overlap]]+1.<0,1>[[Overlap]]";
    Overlap * o=new Overlap(this,ovrOff);
    idx()->setOverlap(o);

    // get the inverse overlap oparator (block-inverse)
    idx()->setInverseOverlap(new Sinverse(this,o));
}



DiscretizationHybrid::Overlap::Overlap(const DiscretizationHybrid *H, string OffDiag):
    OperatorTree("Overlap(Hybrid)",H->idx(),H->idx())
{
    _block.resize(H->comp.size());
    for(unsigned int k=0;k<H->comp.size();k++){
        _block[k].resize(H->comp.size(),0);
        for(unsigned int l=0;l<H->comp.size();l++){
            if(k==l){
                const OperatorTree* ovTree=dynamic_cast<const OperatorTree*>(H->comp[k]->idx()->overlap());
                if(ovTree==0)DEVABORT("overlap is not OperatorTree: "+H->comp[k]->idx()->strNode());
                _block[k][l]=ovTree;
            }
            else{
                OperatorDefinition defkl(OperatorData::extractBlock(OffDiag,k,l),"");
                UseMatrix mult=UseMatrix::Constant(2,2,1.);
                defkl=defkl.constrain(mult,H->idx(),H->idx());
                if(defkl.str()!=""){
                    OperatorTree* ot = new OperatorTree(defkl.str(),defkl.str(),H->comp[k]->idx(),H->comp[l]->idx());
                    ParallelOperator::bcast(ot);
                    _block[k][l]=ot;
                }
            }
            //HACK - this may not create a legitimate operatorTree
            if(_block[k][l]!=0)childAdd(const_cast<OperatorTree*>(_block[k][l]));
        }
    }

    // in the range of the extra basis, all operators will be computed exactly
    // HACK: overlap will NOT be recomputed, i.e. there is an inconsistency between overlap and other operators
    for(const OperatorTree* offS=_block[0][1];offS!=0;offS=offS->nodeNext(_block[0][1])){
        BasisDVR* bd=const_cast<BasisDVR*>(dynamic_cast<const BasisDVR*>(offS->jIndex->basis()));
        if(bd!=0)bd->setDVR(false);
    }

    if(not H->comp[1]->idx()->isOverlapDiagonal()){
        if(not ReadInput::main.flag("DEBUGfem","run FEM basis - compute exact integrals everywhere"))
            DEVABORT("basic overlap not diagonal - implement proper inverse");
    }
}

void DiscretizationHybrid::Overlap::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    Y.scale(B);
    B=1.;
    for(size_t k=0;k<_block.size();k++)
        for(size_t l=0;l<_block[k].size();l++)
            if(_block[k][l]!=0){
                _block[k][l]->apply(A,*Vec.child(l),B,*Y.child(k));
            }
}

DiscretizationHybrid::Sinverse::Sinverse(DiscretizationHybrid *H, const Overlap *Ovr)
    :Inverse("InvOvr:"+H->name,H->idx(),H->idx()),sInvBA(0),sAB(0){

    if(H->comp.size()!=2)ABORT("Sinverse only for exactly two components in Hybrid");
    if(H->comp[0]->idx()->sizeStored()>H->comp[1]->idx()->sizeStored())
        PrintOutput::warning("sort summands such that first is smaller");

    sAinv=H->comp[0]->idx()->inverseOverlap();
    sBinv=H->comp[1]->idx()->inverseOverlap();

    // block maps
    sAB=Ovr->block(0,1);

    aVec=new Coefficients(H->comp[0]->idx());
    bVec=new Coefficients(H->comp[1]->idx());
    aVec->treeOrderStorage();
    bVec->treeOrderStorage();

    zInv=Eigen::MatrixXcd::Zero(aVec->size(),aVec->size());
    if(Ovr->block(0,1)==0 and Ovr->block(1,0)==0){
        PrintOutput::warning("hybrid blocks seem to be orthogonal");
    }

    Coefficients aJ(*aVec);
    Coefficients bTmp(*bVec);
    aJ.treeOrderStorage();
    if(sAB!=0)sInvBA=new OperatorVectors("SbInvBA",Ovr->block(1,0)->iIndex,Ovr->block(1,0)->jIndex);
    for(int j=0;j<aVec->idx()->sizeCompute();j++){
        aJ.setToZero();
        aJ.data()[j]=1.;// j'th unit vector ej
        Ovr->block(0,0)->apply(1.,aJ,0.,*aVec);  // a = Sa ej
        if(sAB!=0){
            Ovr->block(1,0)->apply(1.,aJ,0.,*bVec); // b = C ej
            sBinv->apply(1.,*bVec,0.,bTmp);
            sInvBA->insertColumn(j,bTmp);
            sAB->apply(-1.,bTmp,1.,*aVec); // a <- a - C^H Sb^-1 C ej = [Sa - C^H Sb^-1 C] ej
        }
        zInv.col(j)=Eigen::Map<Eigen::MatrixXcd>(aVec->data(),zInv.rows(),1);
    }
    zInv=zInv.inverse();
}
void DiscretizationHybrid::Sinverse::apply(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{
    //if(MPIwrapper::Size()>1)ABORT("for now, only single processor");

    Coefficients* yPtr;
    const Coefficients* vecPtr;

    if(MPIwrapper::Size() > 1) {
        Parallel::allGather(gX.get(), const_cast<CoefficientsLocal*>(dynamic_cast<const CoefficientsLocal*>(&Vec)));
        vecPtr = gX.get();
        yPtr = gY.get();
    }
    else {
        yPtr = &Y;
        vecPtr = &Vec;
    }
    apply(A,*dynamic_cast<const Coefficients*>(vecPtr),B,*dynamic_cast<Coefficients*>(yPtr));
    if(MPIwrapper::Size() > 1)
        Parallel::scatter(gY.get(), dynamic_cast<CoefficientsLocal*>(&Y), MPIwrapper::master());
}

void DiscretizationHybrid::Sinverse::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    // use Woodbury formula, see tsurff.pdf for the algorithm
    if(Y.childSize()!=2)ABORT("hybrid basis must have exactly 2 components");

    Y.scale(B); // res = B*Y;
    if(A!=1.)ABORT("not for A!=1");

    *aVec=*Vec.child(0); // a = va
    sBinv->apply(1.,*Vec.child(1),0.,*bVec); // b = Sb^-1 vb
    if(sAB!=0)sAB->apply(-1.,*bVec,1.,*aVec); // a = va - C^T b
    Eigen::Map<Eigen::VectorXcd>(aVec->data(),zInv.rows())=zInv*Eigen::Map<Eigen::VectorXcd>(aVec->data(),zInv.rows());// a <- Z^-1 a
    if(sAB!=0)sInvBA->apply(-1.,*aVec,1.,*Y.child(1));  // res=res - Sb^-1 C a = res - Sb^-1 C Z^-1 (va - C^T Sb^-1 vb)

    Y.child(0)->operator+=(*aVec); // res=res + a = res + Z^-1(va-C^T Sb^-1 b)
    Y.child(1)->operator+=(*bVec); // res=res + b = res + Sb^-1 vb
}

void DiscretizationHybrid::Sinverse::parallelSetup() const {
    //if(MPIwrapper::Size()>1)ABORT("for now, only single processor");
    if(MPIwrapper::Size()>1) {//ABORT("for now, only single processor");
        gX.reset(new CoefficientsGlobal(iIndex));
        gY.reset(new CoefficientsGlobal(jIndex));
        dynamic_cast<const Inverse*>(sAinv)->parallelSetup();
        dynamic_cast<const Inverse*>(sBinv)->parallelSetup();
    }
}
