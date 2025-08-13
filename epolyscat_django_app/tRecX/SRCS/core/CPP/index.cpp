// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "index.h"

#include "memInfo.h"
#include <algorithm>
#include "basisOrbital.h"
#include "indexNew.h"
#include "indexAssembled.h"
#include "indexExtract.h"
#include "discretization.h"
#include "discretizationDerived.h"
#include "coefficients.h"
#include "coefficientsFloor.h"
#include "parameters.h"
#include "operatorTensor.h"
#include "readInput.h"
#include "printOutput.h"
#include "overlapDVR.h"
#include "inverseDvr.h"
#include "inverseDVRmixFE.h"
#include "inverseFem.h"
#include "operatorTree.h"
#include "operatorFloor.h"
#include "str.h"
#include "basisNdim.h"
#include "operatorDefinition.h"
#ifdef _USE_HACC_
#include "haccInverse.h"
#endif
#include "tRecXchecks.h"
#include "basisMat1D.h"
#include "basisMat2D.h"
#include "basisSub.h"
#include "parallelOperator.h"
#include "basisDvr.h"
#include "basisGrid.h"

#include "indexConstraint.h"
#include "indexOverlap.h"

#include "coefficients.h"
#include "parallel.h"

#include "basisVector.h"

#include "blockView.h"
#include "asciiFile.h"


using namespace std;
using namespace tools;

int cnt=0;

bool Index::noDum=true;



unsigned int Index::npos=INT_MAX; //!< indicate undefined index counter

static std::map<std::string,size_t> memoryUsage;
std::vector<std::string> Index::_allAxisNames(1,"NONE");
std::map<std::string,std::string> Index::_axisSubset;

std::vector<const BasisAbstract*> Index::_allBases(1,0);
const BasisAbstract* Index::basis() const {return _allBases[_indexBas];}
void Index::setBasis(const BasisAbstract* Bas){
    vector<const BasisAbstract*>::const_iterator b=std::find(_allBases.begin(),_allBases.end(),Bas);
    if(b==_allBases.end()){
        if(_allBases.size()>UINT16_MAX)ABORT("extremely large number of different bases (option -D_LARGE_ to be introduced");
        _indexBas=_allBases.size();
        _allBases.push_back(Bas);
    } else {
        _indexBas=b-_allBases.begin();
    }
}

OperatorTree* findNodeAt(OperatorTree* opNode, const std::vector<unsigned int>& path) {
    // jump over levels created by the summands in more complex overlap operators
    // these nodes do not change the associated index
    // opNode->nodeAt(path) merely follows the path amongst its children
    // HACKY
    // TODO use a more robust code
    
    if(opNode->childSize()==0) {
        if(path.size()==0) return opNode;
        else return 0;
    } else
        if(opNode->child(0)->iIndex==opNode->iIndex) {//skip
            if(opNode->childSize() != 1)DEVABORT("can't skip branching tree level");
            return findNodeAt(opNode->child(0), path);
        } else { // descend
            if(path.size() == 0) return opNode;
            if(path[0] >= opNode->childSize()) return 0;
            return findNodeAt(opNode->child(path[0]), std::vector<unsigned int>(path.begin()+1,path.end()));
        }
}

std::set <const Index*> _hybridList;

long long axisNameCount=0;
bool Index::build=false;
const std::string Index::axisName() const {
    return _allAxisNames[_indexAx];
}
const std::string Index::getAxisSubset(std::string AxisName){
    auto pSub=_axisSubset.find(AxisName);
    if(pSub==_axisSubset.end())return "";
    return pSub->second;
}

const std::string Index::axisSubset() const {
    //    auto pSub=_axisSubset.find(axisName());
    //    if(pSub==_axisSubset.end())return "";
    //    return pSub->second;
    return getAxisSubset(axisName());
}
void Index::setAxisSubset(string Subset){
    if(axisName()=="NONE")DEVABORT("set axis name first");
    // register non-empty Subset (avoid erase for element axes)
    if(Subset!="")_axisSubset[axisName()]=Subset;
}

void Index::setAxisName(string Name){
    std::vector<std::string>::const_iterator a=find(_allAxisNames.begin(),_allAxisNames.end(),Name);
    if(a==_allAxisNames.end()){
        if(_allAxisNames.size()>UINT8_MAX)DEVABORT(Str("more than")+int(UINT8_MAX)+"different axisName's (implement -D_LARGE_");
        _indexAx=_allAxisNames.size();
        _allAxisNames.push_back(Name);
    }
    else
        _indexAx=a-_allAxisNames.begin();
}


OperatorTree *Index::localOverlap() const {return IndexOverlap::getTree(this,0);}
OperatorTree * Index::localInvOvr() const {return IndexOverlap::getTree(this,1);}

static std::map<const Index*,const Inverse*> _inverseOverlapList;
void Index::setInverseOverlap(const Inverse *Inv){
    IndexOverlap::setInverseOverlap(this,Inv);
}

const Inverse * Index::inverseOverlap() const {
    if(not IndexOverlap::inverseOverlap(this)){
        PrintOutput::DEVwarning(Sstr+"missing inverse overlap for "+hierarchy()+this+" - try build now (should fix constructor)");
        const_cast<IndexNew*>(dynamic_cast<const IndexNew*>(this))->buildOverlap();
    }
    return IndexOverlap::inverseOverlap(this);
}


Index::Index(const Index &Other){
    nodeCopy(&Other,false);
    for(size_t k=0;k<Other.childSize();k++)
        childAdd(new Index(*Other.child(k)));
    sizeCompute();
}

Index::Index():_size(npos),_indexBas(0),_indexAx(0),_indexKind(' '){}
Index::Index(const std::vector<const BasisAbstract *> Bases, const std::vector<std::string> Names, unsigned int FloorLevel)
    :Index()
{
    setAxisName("NONE");
    if(Bases.size()==0){
//        setBasis(BasisSet::getDummy(1));
        setBasis(BasisVector::factory("Vector:1"));
    }
    else{
        setBasis(Bases[0]);
        setAxisName(Names[0]);
        for(size_t k=0;k<Bases[0]->size();k++)
            childAdd(new Index(vector<const BasisAbstract*>(Bases.begin()+1,Bases.end()),
                               vector<string>(              Names.begin()+1,Names.end())));
        if(parent()==0){
            if(FloorLevel>Bases.size())
                resetFloor(Bases.size()-1);
            else
                resetFloor(FloorLevel);
        }
    }
    sizeCompute();
}

Index::~Index() {
    // cleanup the registries

}

//addNonEmpty(this,Constraint,Path,Pos,removed);
static void addNonEmpty(Index* Idx, const vector<Axis> & Ax, const IndexConstraint* Constraint,std::vector<unsigned int>  Pos,
                        std::vector<const Index *> Path,std::vector<int> & removed)
{
    if(Constraint !=0 and not Constraint->includes(Path, Pos)){
        removed.push_back(Pos.back());
    }
    else {
        Index* child(0);
        try{child=new Index(Ax,Constraint,Pos,Path);
        }catch(Index::empty_subtree_exception& ex){
            delete child;
            child=0;
        }
        if(child != 0){
            Idx->childAdd(child);
        }else{
            removed.push_back(Pos.back());
            delete child;
        }
    }

}

static void setBasisRemoved(Index* Idx,const std::vector<int> & Removed){
    if(Removed.size() == Idx->basis()->size()){
        throw Index::empty_subtree_exception();
    }
    if(Removed.size() > 0){
        if(dynamic_cast<const BasisDVR*>(Idx->basis()))ABORT("cannot constrain DVR basis (for now)");
        Idx->setBasis(Idx->basis()->remove(Removed));
    }
}

Index::Index(const vector<Axis> & Ax, const IndexConstraint* Constraint, std::vector<unsigned int> Pos, std::vector<const Index *> Path)
    :Index(){

    if(Ax.size()>Path.size())
    {
        Pos.push_back(0);
        Path.push_back(this);

        setAxisName(Ax[Path.size()-1].name);

        if(Ax[Path.size()-1].basDef.size()==1 or Ax[Path.size()-1].bases.size()==1){

            if(Ax[Path.size()-1].bases.size()==1){
                setBasis(Ax[Path.size()-1].bases[0]);
            } else {
                BasisSetDef curDef=Ax[Path.size()-1].basDef[0].resolveDependence(Pos,Path);
                IndexNew::resolveDependence(curDef,Pos,Path);
                setBasis(BasisAbstract::factory(curDef));
            }

            if(basis()->size()==0)ABORT("zero-size basis on "+strNode());

            std::vector<int> removed;
            for(unsigned int k=0;k<basis()->size();k++){
                Pos.back()=k;
                addNonEmpty(this,Ax,Constraint,Pos,Path,removed);
            }
            setBasisRemoved(this,removed);
        }

        else if(Ax[Path.size()-1].basDef.size()>1){
            // finite element axis
            setBasis(BasisAbstract::factory("Vector:"+tools::str(Ax[Path.size()-1].basDef.size())));
            vector<Axis> ax(Ax);

            const Axis* fA=&Ax[Path.size()-1];
            std::vector<int> removed;
            for(unsigned int k=0;k<basis()->size();k++){
                // append single element floor axis
                Pos.back()=k;
                BasisSetDef bDef=fA->basDef[k];
                if(Constraint!=0){
                    // remove first/last FE function, if needed (very clumsy)
                    if(Pos.back()>0){
                        Pos.back()--;
                        if(not Constraint->includes(Path,Pos))bDef.first=true;
                        Pos.back()++;
                    }
                    if(Pos.back()<basis()->size()-1){
                        Pos.back()++;
                        if(not Constraint->includes(Path,Pos))bDef.last=true;
                        Pos.back()--;
                    }
                }
                ax.push_back(Axis(fA->name,fA->comsca,bDef));
                if(Constraint!=0 and not Constraint->includes(Path,Pos)){
                    // skip if outside constraint
                    removed.push_back(k);
                    DEVABORT("constraint on FE axis not working yet");
                } else {
                    childAdd(new Index(ax,Constraint,Pos,Path));
                }
                ax.pop_back();
            }
            if(removed.size() == basis()->size()){
                throw empty_subtree_exception();
            }
            if(removed.size() > 0){
                setBasis(basis()->remove(removed));
            }

        }
        Pos.pop_back();
        Path.pop_back();
    }

    else
    {
        if(Ax.size()<Path.size())DEVABORT("this must not happen");
        // end of axis hierarchy
        setBasis(BasisAbstract::factory("Vector:1"));
        setAxisName("NONE");
        unsetFloor();
    }

    if(Path.size()==0){
        sizeCompute();

        //HACK bad hacking for floor positions: looks like inconsistent use of floor concept in hacc
        if(Ax.size()==2 and Ax[0].coor.name()=="Vec" and Ax[1].coor.name()=="Ion"){
            resetFloor(1);
        }
        else if(Ax.back().coor.name()=="Ion" or Ax.back().coor.name()=="Neutral")
            resetFloor(0); //haCC uses full basis
        else if(descend(Ax.size())->isLeaf())
            resetFloor(Ax.size()-1); // no finite elements, take last axis as floor
        else
            resetFloor(Ax.size());

        // compute overlap matrices in places where it makes sense
        IndexOverlap::localOverlapAndInverse(this,0,0);

        if(isOverlapDiagonal()){
            setOverlap(new OverlapDVR(this));
            setInverseOverlap(Inverse::factory(this));
        }
        else
            setInverseOverlap(Inverse::factory(this));
        DEVABORT("obsolet for "+strNode());

    }
}

const BasisAbstract * Index::basisFromGrid(int ContractFactor, std::vector<const Index*> Path) const {
    const BasisGrid *g=basis()->grid();
    if(g==0)return basis();

    Coordinate coor=Coordinate::fromString(axisName());
    if(coor.defaultFunction()=="useIndex")DEVABORT("\"useIndex\" is obsolete, replace by \"vector\"");
    if(coor.defaultFunction()=="vector")return basis();
    if(not subEquivalent())return basis();

    Str db(coor.defaultFunction(),"");
    vector<double>par;
    if(tools::findFirstOutsideBrackets(db,"{","[","]")!=string::npos){
        string depAx=tools::stringInBetween(db,"{","}"); // dependent axis
        if(axisName().find_first_of("0123456789")!=string::npos)
            depAx+=axisName().substr(axisName().find_first_of("0123456789")); // append numbering (if any)
        for(const Index* depIdx: Path){
            if(depIdx->axisName()==depAx){
                db=Str(db.substr(0,db.find("{")),"")+"["+depIdx->basis()->physical(depIdx->childSize())+"]";
                par.push_back(depIdx->basis()->physical(depIdx->childSize()));
            }
        }
    }
    int order=basis()->size()/ContractFactor;
    BasisSetDef def(order,coor.qmin,coor.qmax-coor.qmin,coor.defaultFunction(),true,true,true,
                    coor,{0,order-1},ComplexScaling(),false,par);
    return BasisAbstract::factory(def);
}


Index* Index::toIndexBasis(std::vector<int> ContractFactor,std::vector<const Index*> Path) const
{
    Index* idx=new Index();
    idx->nodeCopy(this,false);

    if(isLeaf())return idx;

    if(ContractFactor.size()<=Path.size())ContractFactor.resize(Path.size()+1,1);
    idx->setBasis(basisFromGrid(ContractFactor[Path.size()],Path));
    Path.push_back(idx);
    for(size_t k=0;k<idx->basis()->size();k++)
        idx->childAdd(child(k)->toIndexBasis(ContractFactor,Path));

    if(Path.size()==1)idx->sizeCompute();
    return idx;
}

bool Index::isOverlapDiagonal() const {
    if(sizeStored()==1)return true; // trivially diagonal
    if(localOverlap()==0)ABORT("set up s0 before testing for it's being diagonal");
    if(localOverlap()->iIndex!=localOverlap()->jIndex){
        Sstr+localOverlap()->iIndex+localOverlap()->jIndex+Sendl;
        Sstr+localOverlap()->iIndex->root()+localOverlap()->jIndex->root()+Sendl;
        Sstr+localOverlap()->iIndex->root()+str()+Sendl;
        Sstr+localOverlap()->iIndex->strNode()+Sendl;
        Sstr+localOverlap()->jIndex->strNode()+Sendl;
        DEVABORT("this should never happen");
    }
    if(hasFloor()){
        if(not localOverlap()->isLeaf())return localOverlap()->descend()->floor()->isDiagonal();
        //HACK if dummy, assume set up on different process
        if(localOverlap()->str().find("DUM")!=string::npos)return true;
        return localOverlap()->floor()->isDiagonal();
    }
    for(size_t k=0;k<childSize();k++){
        if(not child(k)->isOverlapDiagonal())
            return false;
    }
    return true;
}

void Index::testInverseOverlap() const {

    if(tRecX::off("inverseOverlap")
            or ReadInput::main.flag("DEBUGfem","run FEM basis - compute exact integrals everywhere")
            )return;

    // Note: test is per coefficient, this number should be stringent
    double eps=1e-12;

    Coefficients psi1(this);
    Coefficients psi2(this);
    Coefficients psi3(this);
    psi1.setToRandom();

    // make continuous and pseudo-normalize
    psi1.makeContinuous();
    psi1*=1.0/psi1.norm();

    for(int k=0;k<2;k++){
        psi1.makeContinuous();
        psi1.scale(1./psi1.norm());
        overlap()->apply(1.,psi1,0.,psi2);
        inverseOverlap()->apply(1.,psi2,0.,psi3);

        if(k<1)psi1=psi3; // first run acts as projection
        psi3-=psi1;
    }
    double error1=psi3.norm();

    // Sinv S psi and S Sinv psi
    for(int k=0;k<2;k++){
        psi1.makeContinuous();
        psi1.scale(1./psi1.norm());
        inverseOverlap()->apply(1.,psi1,0.,psi2);
        overlap()->apply(1.,psi2,0.,psi3);
        psi3.makeContinuous();
        if(k==0)psi1=psi3; // first run acts as projection
        psi3-=psi1;
    }
    double error2 = psi3.norm();

    double tol = isHybrid()? 1.e4 : 1.e3; // hybrid needs a higher tolerance
    Str mess("","");
    mess=mess+" inverse overlap: S Sinv, error="+tools::str(error1,2)+" ... Sinv S, error="+tools::str(error2,2)+" ("+hierarchy()+")";
    if (max(error1,error2)> eps){
        PrintOutput::DEVwarning(Str("Error")+mess);
        if(error1+error2>eps*tol){
            mess=Sstr+MPIwrapper::Rank()+mess+"\nInverse overlap error outside tolerance\n"+psi3.str();
            if(isHybrid())mess+"\n check basis for over-completeness";
            ABORT(mess);
        }
    } else
        PrintOutput::DEVmessage(Str("OK")+mess);

    if(not overlap()->isComplexSymmetric(1e-12))DEVABORT("not symmetric");

}

const OperatorAbstract* Index::ovrNonSingular() const {
#ifdef _USE_HACC_
    const HaccInverse* inv=dynamic_cast<const HaccInverse*>(inverseOverlap());
    if(inv==0)return overlap();
    return inv->ovrNonSingular();
#else
    return overlap();
#endif
}



// reconstruct Index path to present Index
vector<const Index*> Index::path() const {
    if(parent()==0)return vector<const Index*>(0);
    vector<const Index*> p(parent()->path());
    p.push_back(parent());
    return p;
}

unsigned int Index::heightAboveFloor() const{
    if (hasFloor()){ return 0; }
    else {return child(0)->heightAboveFloor()+1; }
}
unsigned int Index::heightAboveBottom() const{
    if (isBottom() or isLeaf()) {return 0;} // some Index's may not have proper Bottom
    else { return child(0)->heightAboveBottom()+1; }
}


void Index::boundaryIndices(std::vector<unsigned int> &Boundaries) const{
    vector<unsigned int> iglob(contractedNumbering());

    // duplicate global indices are boundary indices
    Boundaries.clear();
    for(unsigned int k=0;k<globalLength();k++){
        vector<unsigned int>::iterator it0=find(iglob.begin(),iglob.end(),k);
        if(it0==iglob.end())ABORT("index missing from iglob");
        if(find(it0+1,iglob.end(),k)!=iglob.end())Boundaries.push_back(k);
    }
}

unsigned int Index::globalLength() const{
    vector<unsigned int> iglob(contractedNumbering());
    return *std::max_element(iglob.begin(),iglob.end())-*std::min_element(iglob.begin(),iglob.end())+1;
}

void Index::globalElementBoundary(double ElemBound, std::vector<unsigned int> &GlobalIndex) const {

    const Index* iFem=this;
    while(iFem!=0 and iFem->continuity()==Index::npos)iFem=iFem->descend();
    if(iFem==0)ABORT("there is no continuity axes - cannot get inter-element boundary");

    // map from index position to global matrix position
    vector<unsigned int> iglob(contractedNumbering());
    for(const Index* iCont=iFem->descend(iFem->continuity()-iFem->depth());iCont!=0;iCont=iCont->nodeRight(this))
        if(abs(iCont->basis()->integrable()->lowBound()-ElemBound)<1.e-14)
            GlobalIndex.push_back(iglob[iCont->posIndex(this)]);

    if(GlobalIndex.size()==0)ABORT(Str("could not locate coordinate")+ElemBound+"on axis"+iFem->axisName());

    // exclude further continuity levels
    while(0!=(iFem=iFem->descend()) and iFem->continuity()==Index::npos);
    if(iFem!=0)ABORT("not for multiple continuity levels");
}

void Index::matrixContract(const UseMatrix &Mat,UseMatrix & GMat) const {
    // get global indices
    vector<unsigned int>iglob(contractedNumbering());

    // check matrix dimensions
    if(iglob.size()!=Mat.rows() or iglob.size()!=Mat.cols())
        ABORT("input matrix dimensions do not match index size "
              +tools::str(int(iglob.size()))+" X "+tools::str(int(iglob.size()))+" vs. "+tools::str(Mat.rows())+" "+tools::str(Mat.cols()));

    // contract columns
    UseMatrix cmat=UseMatrix::Zero(Mat.rows(),globalLength());
    for (unsigned int j=0;j<Mat.cols();j++)
        cmat.block(0,iglob[j],Mat.rows(),1)+=Mat.col(j);

    // contract rows
    GMat=UseMatrix::Zero(globalLength(),globalLength());
    for (unsigned int i=0;i<cmat.rows();i++)
        GMat.block(iglob[i],0,1,cmat.cols())+=cmat.row(i);
}
vector<complex<double> > Index::dvrWeigContract() const {

    vector<complex<double> > diag;
    const Index * l=firstFloor();
    for(;l!=0;l=l->nodeRight()){
        vector<complex<double> > ovr;
        l->dvrWeights(ovr);
        diag.insert(diag.end(),ovr.begin(),ovr.end());
    }
    vector<unsigned int> iglob(contractedNumbering());
    vector<complex<double> > ret(globalLength(),0.);
    for(unsigned int i=0;i<iglob.size();i++)ret[iglob[i]]+=diag[i];
    return ret;
}


std::vector<double> Index::grid(string Axis) const {
    // descend to Axis
    vector<double> g;
    const Index * iAx=this;
    while(iAx!=0 and iAx->axisName()!=Axis)iAx=iAx->descend();
    if(iAx==0)return g;

    if(iAx->continuity()!=Index::npos)iAx=descend(iAx->continuity());
    if(not iAx->basis()->grid())return g;

    // run through all grid elements
    while(iAx!=0){
        for(size_t k=0;k<iAx->basis()->grid()->size();k++){
            g.push_back(iAx->basis()->grid()->mesh()[k]);
        }
        iAx=iAx->rightSibling();
    }
    return g;
}


// user interfaces for contracted numbering
vector<unsigned int> Index::contractedNumbering() const {
    vector<unsigned int> numb,mult;
    int curMax=-1,curPos=-1;
    contractedNumbering(numb,mult,this,curPos,curMax);
    return numb;
}

void Index::multiplicities(std::vector<unsigned int> &Glob, std::vector<double>& Norm) const{
    vector<unsigned int> mult;
    contractedNumbering(Glob,mult);
    for(size_t k=0;k<mult.size();k++)Norm.push_back(1./sqrt(double(mult[k])));
}


void Index::contractedNumbering(std::vector<unsigned int> & Numb,std::vector<unsigned int> & Mult) const{
    int cur=-1,pos=-1;
    contractedNumbering(Numb,Mult,this,cur,pos);
}

void Index::bottomExpandAll() const {
    if(not firstLeaf()->isBottom())return;
    for(Index * ix=firstLeaf();ix!=0 and ix->isBottom();ix=ix->nodeRight())ix->bottomExpand();
}

void Index::bottomUnexpandAll() const {
    Index* ix=firstLeaf();
    while(ix and not ix->isBottom())ix=const_cast<Index*>(ix->parent());
    for(;ix and ix->isBottom();ix=ix->nodeRight())ix->bottomUnexpand();
}


void Index::bottomExpand() const {
    if(not isBottom() or not isLeaf())return;
    for(size_t k=childSize();k<basis()->size();k++){
        const_cast<Index*>(this)->leafAdd();
        childBack()->_indexKind='X';
        childBack()->_size=1;
    }
}
void Index::bottomUnexpand() const {
    if(not isBottom() or
            not childSize() or
            child(0)->_indexKind!='X')
        return;
    for(int k=childSize()-1;k>=0;k--)
        const_cast<Index*>(this)->childErase(k);
}


// internal interface
void Index::contractedNumbering(std::vector<unsigned int> & Numb,std::vector<unsigned int> & Mult, const Index* From,int&CurMax,int &CurPos) const{

    if(this==From){
        Numb.assign(From->sizeStored(),From->sizeStored());
        Mult.assign(Numb.size(),1);
    }

    bottomExpand();
    if(isLeaf()){
        // assign, unless it has been assigned previously
        CurPos++;
        if(Numb[CurPos]==From->sizeStored())Numb[CurPos]=++CurMax;
    }else{
        if(basis()->size()!=childSize()){
            DEVABORT("not for subtrees");
        }

        const Index* lN=lowerNeighbor();
        if(lN)lN->bottomExpand();
        for(unsigned int k=0;k<childSize();k++){
            // add new "indices"
            child(k)->contractedNumbering(Numb,Mult,From,CurMax,CurPos);
            if(lN!=0 and k==basis()->integrable()->lowerMargin()){
                // has lower neighbor - reset lower margin indices to lower neighbor's upper margin indices
                unsigned int cur=child(k)->posIndex(From);
                unsigned int nei=lN->child(lN->basis()->integrable()->upperMargin())->posIndex(From);
                for(size_t i=0;i<child(k)->sizeStored();i++){
                    Mult[cur+i]*=2;
                    Mult[nei+i]*=2;
                    if(cur>nei)Numb[cur+i]=Numb[nei+i];
                    else       Numb[nei+i]=Numb[cur+i];
                }
                // last index is no longer the largest,re-determine
                if(cur>nei)CurMax=*std::max_element(Numb.begin()+nei,Numb.begin()+cur+child(k)->sizeStored());
            }
        }
        if(lN)lN->bottomUnexpand();
    }
    bottomUnexpand();
}

//void Index::unGlobal(UseMatrix &Eigen,bool AsDual, vector<Coefficients*> & Evec, vector<int> cols) const{
//    vector<unsigned int>cI(contractedNumbering());
//    for (unsigned int col=0; col<cols.size(); col++) {
//        Evec.push_back(new Coefficients(this));
//        const Coefficients * leaf=Evec.back()->firstLeaf();
//        while(leaf!=0){
//            unsigned int pos=leaf->idx()->posIndex(this);
//            for(unsigned int k=0;k<leaf->size();k++){
//                leaf->floorData()[k]=Eigen(cI[pos+k],cols[col]).complex();
//            }
//            leaf=leaf->nextLeaf();
//        }
//        // contracted dual vectors carry the contraced overlap matrix
//        // remove that extra factor 2 from margins
//        Evec.back()->makeContinuous(sqrt(0.5));
//    }
//}

void Index::cleanFloors(const bool FloorFound){
    if(FloorFound)unsetFloor();
    for(size_t k=0;k<childSize();k++)child(k)->cleanFloors(FloorFound or hasFloor());
}

void Index::resetFloor(unsigned int Level){
    if(Level<depth()){
        unsetFloor();
    }
    else if (Level>depth()){
        unsetFloor();
    } else {
        setFloor();
    }
    for(unsigned int k=0;k<childSize();k++)child(k)->resetFloor(Level);
    if(parent()==0)sizeCompute();
}

void Index::setFloorAuto(std::vector<std::string> & FemAxes) {
    // default floor level:
    // if no FE axes are found, bottom-level is floor
    // if FE, highest FE sub-level
    unsetFloor();
    if(isBottom() or std::find(FemAxes.begin(),FemAxes.end(),axisName())!=FemAxes.end())
        setFloor();
    else {
        bool fem=isFem();
        if(fem)FemAxes.push_back(axisName());
        for(size_t k=0;k<childSize();k++)child(k)->setFloorAuto(FemAxes);
        if(fem)FemAxes.pop_back();
    }
}

unsigned int Index::depthInFloor() const
{
    int floorLevel = 0;
    for(const Index* s=this; not s->hasFloor(); s=s->parent()){
        if(s->isRoot()) return npos;
        floorLevel++;
    }
    return floorLevel;
}

void Index::setFloor(unsigned int Level){

    if(continuity()!=Index::npos and Level<depth())ABORT("cannot place floor at or above any finite element level");

    // check for previous occurance of floor
    const Index * f=this;
    while(f!=0){
        if(f->hasFloor())ABORT("found previous floor on "+tools::str(f->index()));
        f=f->descend();
    }

    if(depth()<Level){
        if(hasFloor())ABORT("floor had been defined previously");
        for(unsigned int k=0;k<childSize();k++)child(k)->setFloor(Level);
    }
    else {
        setFloor();
    }
}

bool Index::isBottom() const {
    if(isLeaf())return axisName()!="NONE";
    return const_cast<Index*>(this)->childRef(0)->axisName()=="NONE";
}

bool Index::hasFloor() const{
    return _indexKind=='F';
} //!< is floor index

void Index::setFloor(){
    _indexKind='F';
}

void Index::unsetFloor(){if(_indexKind=='F')_indexKind=' ';}
void Index::setKind(const char Kind){_indexKind=Kind;}

Index* Index::axisIndex(const string Name) const {
    if(axisName()==Name)return const_cast<Index*>(this);
    if(isHybrid()){
        for(size_t k=0;k<childSize();k++){
            Index* idx=child(k)->axisIndex(Name);
            if(idx)return idx;
        }
    }
    if(isLeaf())return 0;
    return descend()->axisIndex(Name);
}

void Index::leafAdd(){
    childAdd(new Index());
    childBack()->setBasis(BasisAbstract::factory("Vector:1"));
    childBack()->setAxisName("NONE");
    childBack()->_size=1;
}


//static string failureCompatible; // file-wide variable (not pretty)
bool Index::subEquivalent(string Mess) const {
    failureCompatible="";
    for (unsigned int k=1;k<childSize();k++)
        if(not child(0)->treeEquivalent(child(k))){
            if(Mess!="")cout<<Mess<<": "<<failureCompatible<<endl;
            return false;
        }
    return true;
}

Index* Index::factor(const std::vector<int> Levels, bool Exclude) const{
    Index* idx=Tree::factor(Levels,Exclude);
    idx->sizeCompute();
    return idx;
}

bool Index::nodeEmpty() const{
    return basis()==0 and not hasFloor();
}

void Index::purge(unsigned int Height){
    if(Height==1)return;
    vector<int> removed;
    vector<int> subset;
    for(int k=childSize()-1;k>=0;k--){
        child(k)->purge(Height-1);
        if(child(k)->isLeaf()){
            childErase(k);
            removed.insert(removed.begin(),k);
        }
        else
            subset.insert(subset.begin(),k);
    }
    if(removed.size()>0 and removed.size()!=basis()->size()){
        setBasis(BasisAbstract::factory(BasisSub::strDefinition(basis(), subset)));
    }
}

void Index::nodeCopy(const Index *Node, bool View){
    if(View)DEVABORT("Index cannot be View=true");
    setBasis(Node->basis());
    setAxisName(Node->axisName());
    _size=npos; // as position in tree is undefined, size is undefined
    _indexKind=Node->_indexKind;
}

bool Index::nodeEquivalent(const Index *Other) const{

    // quick and likely failed checks first
    failureCompatible="axes: "+axisName()+" != "+Other->axisName();
    if(axisName()!=Other->axisName())return false;

    if(basis()!=Other->basis()){
        if(basis()==0 or Other->basis()==0)
            ABORT(Str("undefined basis this:")+basis()+", Other: "+Other->basis()+"at axes"+axisName()+Other->axisName());
        failureCompatible+=":\n"+basis()->strDefinition()+"\n!=\n"+Other->basis()->strDefinition();
        if(not (*basis()==*Other->basis()))return false;
    }

    failureCompatible="floorsA";
    if((hasFloor())!=(Other->hasFloor()))return false;

    failureCompatible="";
    return true;
}


size_t Index::diagnoseSizeOfNode() const {
    size_t siz,sizAll=0;
    if(localOverlap()!=0)    {siz=   localOverlap()->diagnoseSizeOf();memoryUsage["s0"]+=siz;    sizAll+=siz;}

    sizAll           +=sizeof *this;
    if(parent()==0)
        PrintOutput::DEVmessage(Str("Index sizes: ")+memoryUsage["s0"]+memoryUsage["invS0"]+memoryUsage["ovr"]+memoryUsage["invOvr"]+":"+sizAll);
    return sizAll;
}

bool Index::treeEquivalent(const Index * Other, const std::vector<string> OnlyAxes) const{
    if(this==Other)return true;

    failureCompatible+=Sstr+"indexSize"+childSize()+Other->childSize()+"\n"+str(-1)+"\noht\n"+Other->str(-1);
    if(childSize()!=Other->childSize()){
        failureCompatible=Sstr+strNode()+"bottom size"+isBottom() + Other->isBottom() + (size()==Other->size());
        if(not (isBottom() and Other->isBottom() and size()==Other->size()))return false;
    }

    failureCompatible+="not node equivalent";
    if((OnlyAxes.size()==0 or std::find(OnlyAxes.begin(),OnlyAxes.end(),axisName())!=OnlyAxes.end())
            and not nodeEquivalent(Other))return false;

    // equivalence of sub-indices
    for(unsigned int k=0;k<std::min(childSize(),Other->childSize());k++){
        failureCompatible+="subIdx: "+child(k)->strNode()+" : "+Other->child(k)->strNode();
        if(not child(k)->treeEquivalent(Other->child(k),OnlyAxes))return false;
    }
    failureCompatible="";
    return true;
}

/**
 * called in Index-Constructor
 * with arguments (0,0)
 * on root Index.
 * Probably BasicDisc as Discretization class
 * Axis should not be Hybrid, and likely not NDim
 */
void Index::localOverlapAndInverse(OperatorTree* Ovr, OperatorTree* Inv){
    IndexOverlap::localOverlapAndInverse(this,Ovr,Inv);
}

const Index* Index::findAxisStarts(string Start) const{
    const Index *idx=this;
    while(idx!=0 and idx->axisName().find(Start)!=0)idx=idx->nodeNext();
    return idx;
}

unsigned int Index::continuity(unsigned int N) const{
    // descend until continuity level with at least N continuity levels above
    const Index * idx=this;
    while(idx!=0 and (idx->continuity()==Index::npos or idx->contDepth()!=N)){
        idx=idx->nodeNext();
    }
    if(idx==0)return Index::npos;
    return idx->depth()-depth();
}

unsigned int Index::continuity() const{
    // depth of duplicate, if continuity type axisName
    if(Coordinate::isDiscrete(axisName()))return npos;
    return depthOfDuplicate();
}

unsigned int Index::depthOfDuplicate() const{
    if(isLeaf())return npos;
    const Index* s = this;
    do {
        s = s->descend();
        if(not s)break;
        if(s->_indexAx==_indexAx or "k"+s->axisName()==axisName()) return s->depth();
    } while (not s->isBottom() and not s->isLeaf());
    return npos;
}

bool Index::isFem() const {
    if(childSize()==1)return false;
    const Index* idx=descend();
    while(idx!=0 and idx->axisName()!=axisName())idx=idx->descend();
    if(idx==0)return false;
    return idx->basis()->integrable()!=0 or idx->basis()->ndim()!=0;
}

bool Index::isAbsorptive() const {
    if(not basis()->integrable())return false;
    if(basis()->isAbsorptive())return true;
    for(unsigned int k=0;k<childSize();k++)
        if(child(k)->isAbsorptive())return true;
    return false;
}

bool Index::isHybrid() const {
    if(hierarchy().find("&")!=std::string::npos and childSize()>1)return true;
    if(axisName() != "Hybrid") return false;
    if(childSize()<2)return false;
    return child(0)->axisName() != child(1)->axisName();
}

void Index::dvrWeights(std::vector<std::complex<double> >& Weights) const{

    if(not hasFloor())DEVABORT("only on floor level");
    DEVABORT("needs re-implementation");
    Weights.clear();

}

string Index::hierarchy(unsigned int HybridPath) const {
    std::string s=isBottom() or not isLeaf()?axisName():"NONE";
    return isLeaf()?s:s+="."+child(std::min(HybridPath,childSize()-1))->hierarchy();
}

std::string Index::hierarchy_no_NONE() const{
    std::string h=hierarchy();
    return h.substr(h.length()-5)==".NONE"?h.substr(0,h.length()-5):h;
}

std::string Index::coordinates(std::string Hierarchy){
    vector<string> iH=tools::splitString(Hierarchy,'.');
    // remove duplicate coordinate names
    string s;
    for(size_t k=0;k<iH.size();k++){
        if(tools::anyElement(vector<string>(iH.begin()+k+1,iH.end()),tools::equal,iH[k]))continue;
        // special case: a specXXX is followed by kXXX does not count as coordinate
        if(iH[k].find("spec")!=string::npos and Hierarchy.find("k"+iH[k].substr(4))!=string::npos)continue;
        if(iH[k].find("surf")!=string::npos and Hierarchy.find("ValDer"+iH[k].substr(4))!=string::npos)continue;
        if(iH.size()>1 and iH[k]=="NONE")continue; // do not attach a trailing "NONE"
        s+="."+iH[k];
    }
    return s.substr(1); // remove leading "."
}


string Index::coordinates(unsigned int HybridPath) const {
    std::string hier(hierarchy(HybridPath));
    if(basis()->orbital())hier=basis()->orbital()->orbital(0)->idx()->hierarchy();
    std::vector<std::string> ax(tools::splitString(coordinates(hier),'.'));
    std::string res=ax[0];
    for(size_t k=1;k<ax.size();k++)
        if(std::find(ax.begin(),ax.begin()+k,ax[k])==ax.begin()+k)
            res+="."+ax[k];
    return res;
}

unsigned int  Index::sizeCompute() const {
    unsigned int sz=0;

    if(isLeaf()){
        basis()?sz=basis()->size():0;
    }
    for (unsigned int k=0;k<childSize();k++){
        int cSize=child(k)->sizeCompute();
        if(sz+cSize>UINT32_MAX)ABORT("extremely large Index tree (compile flag -D_LARGE_ remains to be introduced)");
        sz+=child(k)->sizeCompute();
    }

    // store, if data is owned
    if(not isView())_size=sz;
    return sz;
}

// total size up to beginning of present Index
unsigned int Index::posIndex(const Index* Root) const {
    if(parent()==0 or this==Root)return 0;
    unsigned int n=parent()->posIndex(Root);
    for(const Index* s=parent()->child(0);s!=this;s=s->rightSibling()){
        n+=s->sizeStored();
    }
    return n;
}

unsigned int Index::axisLevel(const Index *ITree) const
{
    std::string namL=axisName();
    bool conL=continuity()!=npos;
    unsigned int pos=0;
    while(ITree->descend(pos)!=0){
        if(ITree->descend(pos)->axisName()==namL
                and conL==(ITree->descend(pos)->continuity()!=npos))return pos;
        pos++;
    }
    return npos;
}

// total size from beginning of Floor to beginning or present Index
unsigned int Index::posInFloor() const {
    if(hasFloor())return 0;
    if(parent()==0)return Index::npos;
    unsigned int n=parent()->posInFloor();
    if(n<Index::npos){
        for(const Index* s=parent()->child(0);s!=this;s=s->rightSibling()){
            unsigned int sz=s->sizeStored();
            if(sz==Index::npos){
                cout<<"pos "<<sz<<" "<<s->sizeCompute()<<" "<<s->sizeStored()<<" "<<tools::str(s->index())<<" view="<<isView()<<endl;
                ABORT("bad");
            }
            n+=s->sizeStored();
        }
    }
    return n;
}

Index* Index::leafAtPos(unsigned int Pos) const{

    // no more children and not exceeded
    if(childSize()==0)return const_cast<Index*>(this);

    size_t k=0;
    for(;Pos>child(k)->sizeStored();k++)Pos-=Pos>child(k)->sizeStored();
    if(k==childSize())ABORT("Pos >= number of leafs ");
    return child(k)->leafAtPos(Pos);
}

unsigned int Index::nSub(const Index * SubIndex) const {
    if(SubIndex->parent()!=this){
        cout<<"SubIndex, parentIndex, this "<<SubIndex<<" "<<SubIndex->parent()<<" "<<this<<endl;
        ABORT("argument is not sub-index of index");
    }
    unsigned int n=0;
    for(;n<childSize();n++)
        if(child(n)==SubIndex)return n;
    ABORT("argument is not in Idx, childSize="+tools::str(childSize()));
}

std::string Index::strAxes() const {
    std::string res;
    std::string prevNam=axisName();
    size_t prevDepth=depth();
    const Index* idx=this;
    unsigned int minChild=INT_MAX,maxChild=0;
    for(;idx!=0;idx=idx->nodeNext()){
        if(idx->depth()==prevDepth and prevNam==idx->axisName()){
            minChild=std::min(minChild,idx->childSize());
            maxChild=std::max(maxChild,idx->childSize());
        } else {
            res+="\n"+prevNam+"\t"+tools::str(minChild);
            if(minChild<maxChild)res+="..."+tools::str(maxChild);
            if(idx->hasFloor())res+="\tF";
            prevDepth=idx->depth();
            prevNam=idx->axisName();
            minChild=INT_MAX;
            maxChild=0;
        }
    }
    //note: \n appears to be encoded by 2 characters
    return res.substr(1);
}


/// superseed Tree::strData with empty or parents (if any)
string Index::strNode(int Level) const{
    if(Level==Tree_ptrsOnly)return Tree::strNode(Level);

    Str s("("+axisName(),"");
    if(axisSubset()!="" and (parent()==0 or parent()->axisSubset()!=axisSubset()))s+=":"+axisSubset();
    s+=")";
    if(hasFloor()){
        if(MPIwrapper::Size()>1){
            if(Parallel::owner(this)==Parallel::all)s=s+"< a>";
            else if(Threads::isThread(this))s=s+"<th>";
            else if(Parallel::owner(this)==Parallel::none)s=s+"<..>";
            else s=s+" <"+SEP("")+Parallel::owner(this)+SEP("")+">";
        }
        s=s+" F"+sizeStored();
    }
    else          s=s+" "+sizeStored();
    if(isBottom() and childSize()>0)s+="x";
    if(basis())s=s+" "+basis()->str();
    else               s=s+" (-no-basis-)";
    s=s+" p:"+posIndex();
    if(Level==Tree_withPtrs)s=s+" "+this;
    return std::move(s);
}

// write number of branches and floor sizes to file
void Index::writeStructure(std::ofstream &stream) const{
    if(depth()==0){
        tools::write(stream,(int) 1); // storage type code (for use for backward compatibility)
    }
    tools::write(stream,(unsigned int) childSize());
    if(childSize()==0){
        if(not hasFloor()){
            tools::write(stream,(unsigned int) 0);               // empty branch
        } else {
            tools::write(stream,(unsigned int) (sizeCompute()) );// non-empty floor
        }
    } else {
        for (unsigned int n=0;n<childSize();n++)child(n)->writeStructure(stream);
    }
}


/// data header:
/// - list of all axis names, format: _indexAx(1Byte),  lengthOfAxisName(1Byte),   axisName(length in Bytes)
/// - list of all bases definitions: _indexBas(2Bytes),lenghtOfDefinition(2Bytes),basisDefinition(length in Bytes)
///
/// then, recursively:
/// - current Index node's _indexAx,_indexBas,_indexKind,basisSize
/// - basisSize=-1 indicates next level is end-of-Index-tree, otherwise basisSize=basis()->size() (redundand)
/// - next nodes's basisSize child's.
void Index::write(ostream &Stream, bool Enter) const{
    if(Enter){
        int currentCode=11; // code to identify current format
        tools::write(Stream,currentCode);

        { // write list of all axis names in index tree
            std::set<uint8_t>setAxes;
            for(const Index* idx=this;idx!=0;idx=idx->nodeNext()){
                if(setAxes.find(idx->_indexAx)==setAxes.end()){
                    setAxes.insert(idx->_indexAx);
                    tools::write(Stream,idx->_indexAx);
                    string sAx=idx->axisName();
                    int    lAxInt=sAx.length();
                    if(currentCode==10){
                        uint8_t lAx=lAxInt;
                        tools::write(Stream,lAx);
                    }
                    else tools::write(Stream,lAxInt);
                    for(auto c: sAx)tools::write(Stream,&c,sizeof(c));
                }
            }
            uint8_t last=UINT8_MAX;
            tools::write(Stream,last);
        }


        { // write list of all bases in index tree
            std::set<uint16_t>setBases;
            for(const Index* idx=this;idx!=0;idx=idx->nodeNext()){
                if(setBases.find(idx->_indexBas)==setBases.end()){
                    setBases.insert(idx->_indexBas);
                    tools::write(Stream,idx->_indexBas);
                    string sBas=idx->basis()->strDefinition();
                    int lenInt=sBas.length();
                    if(currentCode==10){
                        uint16_t len=lenInt;
                        tools::write(Stream,len);
                    }
                    else tools::write(Stream,lenInt);
                    for(auto c: sBas)tools::write(Stream,&c,sizeof(c));
                }
            }
            uint16_t last=UINT16_MAX;
            tools::write(Stream,last);
        }
    }

    // write present Index
    tools::write(Stream,_indexAx);
    tools::write(Stream,_indexBas);
    tools::write(Stream,_indexKind);

    if(isBottom()){
        int basisSize=-1; // number of children: -1 indicates bottom of index
        tools::write(Stream,basisSize);
    } else{
        int basisSize=basis()->size();
        tools::write(Stream,basisSize);
        for(size_t k=0;k<childSize();k++)child(k)->write(Stream,false);
    }
}

/// for data layout on file, see Index::write
Index::Index(istream &Stream, std::map<uint8_t,uint8_t>AxisNumber, std::map<uint16_t,uint16_t>BasisNumber, bool Rewind):Index(){
    if(BasisNumber.size()==0){
        if(Rewind)Stream.seekg(ios_base::beg); // start at beginning of file
        // add all axis names and cross-reference to new numbering
        int code=-1;
        tools::read(Stream,code);
        if(code!=10 and code!=11)
            ABORT(Str("stream does not contain full Index data, code=")+code
                  +" (proc"+MPIwrapper::Rank(MPIwrapper::worldCommunicator())+")");

        uint8_t iNam;
        for(int k=0;k<UINT8_MAX;k++){
            tools::read(Stream,iNam);
            if(iNam==UINT8_MAX)break;

            int lNamInt;
            if(code==10){
                uint8_t lNam;
                tools::read(Stream,lNam);
                lNamInt=lNam;
            }
            else
                tools::read(Stream,lNamInt);

            string axNam(lNamInt,' ');
            for(auto &c: axNam)tools::read(Stream,&c,sizeof(c));
            setAxisName(axNam);
            AxisNumber[iNam]=_indexAx;
        }

        // add all basis and keep cross-reference for new numbering
        uint16_t iBas;
        for(int k=0;k<UINT16_MAX;k++){
            tools::read(Stream,iBas);
            if(iBas==UINT16_MAX)break;

            int lDefInt;
            if(code==10){
                uint16_t lDef;
                tools::read(Stream,lDef);
                lDefInt=lDef;
            }
            else
                tools::read(Stream,lDefInt);

            string basDef(lDefInt,' ');
            for(auto &c: basDef)tools::read(Stream,&c,sizeof(c));
            setBasis(BasisAbstract::factory(basDef));
            BasisNumber[iBas]=_indexBas;
        }
    }

    // read current hierarchy level
    int basisSize;
    uint8_t read8;
    uint16_t read16;
    tools::read(Stream,read8); _indexAx= AxisNumber[read8 ];
    tools::read(Stream,read16);_indexBas=BasisNumber[read16];
    tools::read(Stream,_indexKind);
    tools::read(Stream,basisSize); // read further levels (or set dummy at end of Index)

    if(basisSize==-1){
        basisSize=basis()->size();
    }
    else {
        if(basisSize!=int(basis()->size()))
            DEVABORT(Sstr+"inconsistent construction of index from file:"+axisName()+" "+basis()->str()
                     +" further="
                     +basisSize+"|"+read8+"|"+read16+"|"+_indexKind);
        for(size_t k=0;k<basis()->size();k++)childAdd(new Index(Stream,AxisNumber,BasisNumber,false));
    }
    sizeCompute();
}

// plot all basis functions on a hierarchy level
void Index::axisPlot(std::string File, int Points, double QLow, double QUp) const {
    if(isLeaf()) DEVABORT("Cannot plot leaf index");

    ofstream out;

    for(unsigned int i=0;i<childSize();i++){
        if(child(i)->basis()->integrable()==0) ABORT("child("+std::to_string(i)+": basisIntegrable()=0 - cannot plot axis");
        double plow=max(QLow,child(i)->basis()->integrable()->lowBound());
        double pup=min(QUp,child(i)->basis()->integrable()->upBound());
        if(child(i)->basis()->integrable()->lowBound()>pup or child(i)->basis()->integrable()->upBound()<plow)continue;
        out.open((File+"_"+child(i)->basis()->name()+std::to_string(i)).c_str());
        out<<"# axis plot for basis "+child(i)->basis()->name()<<endl<<"#"<<endl<<"#X"<<" number of function"<<endl;

        // determine the required precision for the x-axis
        int prec = max((int)log10(abs(plow))+1,(int)log10(abs(pup))+1) + (int)log10(Points)+1;

        // maximum number of basis functions on axis
        unsigned int maxcols=child(i)->basis()->integrable()->order();

        // loop through elements on axis
        for(double xi=plow;xi<=pup;xi+=(pup-plow)*(1.-1.e-15)/(max(1,Points-1)))
        {
            out<<setprecision(prec);
            out<<setw(prec+5)<<xi;
            if(child(i)->basis()->integrable()->lowBound()>xi or child(i)->basis()->integrable()->upBound()<xi)continue;
            // get x-points in range
            // evaluate basis functions on this grid
            complex<double>cxi=xi;
            UseMatrix vals;
            vals=child(i)->basis()->integrable()->val(UseMatrix::UseMap(&cxi,1,1),true);
            // write to file
            out<<setprecision(prec);
            for(unsigned int k=0;k<vals.cols();k++)out<<" "<<setw(prec+5)<<vals(0,k).real();
            for(unsigned int k=vals.cols();k<maxcols;k++)out<<" "<<setw(prec+5)<<0.;
            out<<endl;
        }

        out.close(); out.clear();
    }
}

std::string Index::failureCompatible;

/// true if branch structure and floor sizes on file match present index
bool Index::compatibleFile(std::ifstream &stream, int code) const{
    unsigned int _size;
    if(code==0)tools::read(stream,code);
    switch (code){
    case 1:
        tools::read(stream,_size);
        if(_size!=childSize()){
            failureCompatible+=Sstr+"child"+_size+childSize();
            return false;
        }
        if(childSize()==0){ // end of branch
            tools::read(stream,_size);
            if(not hasFloor()){
                failureCompatible+=Sstr+"size"+(_size==0);
                return _size==0;  // empty branch
            } else {
                failureCompatible+=Sstr+"size"+(_size==sizeStored());
                return _size==sizeStored(); // non-empty floor
            }
        } else { // descend into hierarchy
            for (unsigned int n=0;n<childSize();n++)
                if(not child(n)->compatibleFile(stream,code))return false;
            return true;
        }
    default:
        Index inIdx(stream);
        failureCompatible=Sstr+"equiv"+treeEquivalent(&inIdx)+failureCompatible;
        return treeEquivalent(&inIdx);
    }
}

const Index * Index::firstFloor() const {
    if(hasFloor())return this;
    return descend()->firstFloor();
}

unsigned int Index::contDepth() const {
    if(parent()==0)return 0;
    return parent()->contDepth()+(parent()->continuity()!=Index::npos);
}
const Index* Index::lowerNeighbor(unsigned int D) const {
    if(not hasFloor())DEVABORT("apply only on floor level");
    if(D>height())DEVABORT("D exceeds floor dimensions");
    const Index* neig=descend(D)->lowerNeighbor();
    for(size_t k=0;neig!=0 and k<D;k++)neig=neig->parent();
    return neig;
}
const Index* Index::upperNeighbor(unsigned int D) const {
    if(not hasFloor())DEVABORT("apply only on floor level");
    if(D>height())DEVABORT("D exceeds floor dimensions");
    const Index* neig=descend(D)->upperNeighbor();
    for(size_t k=0;neig!=0 and k<D;k++)neig=neig->parent();
    return neig;
}

const Index* Index::lowerNeighbor() const {
    if(Coordinate::isDiscrete(axisName()))return 0;
    // ascend until axis name matches
    const Index* match=parent();
    while(match!=0 and match->axisName()!=axisName())match=match->parent();

    // no match - no neighbor
    if(match==0)return 0;

    vector<unsigned int> idx(index());
    // first sibling does not have lower neighbor
    if(idx[match->depth()]==0){
        if(not basis()->isPeriodic())return 0;
        idx[match->depth()]=match->childSize();
    }

    idx[match->depth()]--;
    return root()->nodeAt(idx);
}

const Index* Index::upperNeighbor() const {
    if(Coordinate::isDiscrete(axisName()))return 0;
    // ascend until axis name matches
    const Index* match=parent();
    while(match!=0 and match->axisName()!=axisName())match=match->parent();

    // no match - no neighbor
    if(match==0)return 0;

    vector<unsigned int> idx(index());
    // last sibling does not have upper neighbor
    if(idx[match->depth()]==match->childSize()-1){
        if(not basis()->isPeriodic())return 0;
        idx[match->depth()]=-1;
    }

    idx[match->depth()]++;
    return root()->nodeAt(idx);
}

double Index::physical() const{
    if(!parent()) return 0.;
    if(!parent()->basis()) return nSibling();
    return parent()->basis()->physical(nSibling());
}

const Index* Index::from() const{
    if(dynamic_cast<const IndexAssembled*>(this))return dynamic_cast<const IndexAssembled*>(this)->from();
    if(dynamic_cast<const IndexExtract*>(this))  return dynamic_cast<const IndexExtract*>(this)->from();
    return this;
}

