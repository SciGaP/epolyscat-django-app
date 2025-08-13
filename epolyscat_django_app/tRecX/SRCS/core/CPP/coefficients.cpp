// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "log.h"

#include "coefficients.h"
#include "coefficientsSparse.h"

#include <climits>
#include <iostream>
#include <fstream>
#include <regex>

// resolve forward declarations
#include "discretization.h"
#include "discretizationGrid.h"
#include "coefficientsFloor.h"
#include "qtEigenDense.h"
#include "index.h"
#include "wavefunction.h"
#include "timer.h"
#include "printOutput.h"
//#include "basisMat.h"
#include "operatorFloor.h"
#include "derivativeFlat.h"
#include "operatorSingle.h"
#include "parallel.h"
#include "mpiWrapper.h"
using namespace std;
#include "eigenNames.h"
#include "operatorDefinition.h"
#include "basisGrid.h"
#include "inverse.h"
#include "basisDvr.h"
#include "quadratureRule.h"
#include "fft.h"
#include "threads.h"
#include "debugInfo.h"
#include "parallelContinuity.h"

#include "timeCritical.h"

#include "discretizationGrid.h"

#include "readInput.h"
#include "asciiFile.h"
#include "plot.h"
#include "algebra.h"

using namespace tools;

// use this for debug switching
std::string Coefficients::DEBUGstatus="";

// local switch for norm type, may be made input-dependent laters
static int normType=0;

std::set<const Coefficients*> Coefficients::_viewFloor;
std::map<const Coefficients*, std::string> Coefficients::_labels;
std::unordered_map<const Coefficients*, std::vector<std::complex<double> >* > Coefficients::_centralStorage;
static std::map<Coefficients*,std::shared_ptr<ParallelContinuity>> _parallelContinuity;

long unsigned int Coefficients::_labelCurrent=0;
std::string Coefficients::hash() const {

    std::string l;
    if((l=_labels[this])==""){
        l="#"+std::to_string(_labelCurrent++);
        _labels[this]=l;
    }
    return l;
}

Coefficients::Coefficients():_cData(0),_cIndex(0){
    if(timeCritical::coefficients)
        PrintOutput::DEVwarning("Coefficients construct during timeCritical",1000);
    parentRef()=0;
}

complex<double>* Coefficients::storageData() const {
    const Coefficients* f=this;
    while(f!=0 and f->storageSize()==0){
        f=f->parent();
    }

    if(f==0)return 0;
    return const_cast<Coefficients*>(f)->nodeStorage()->data();
}

void Coefficients::replaceStorage(std::vector<std::complex<double> > &Stor){
    clearStorage();
    delete _centralStorage[this];
    _centralStorage[this]=new vector<complex<double> >();
    Stor.swap(*nodeStorage());
}

void Coefficients::clearStorage(){
    nodeStorageClear();
    for(size_t k=0;k<childSize();k++)child(k)->clearStorage();
}

bool Coefficients::isContiguous() const {
    return _centralStorage.count(this)==1;
}

/// (CAUTION: data-pointers will be re-directed to new storage)
void Coefficients::treeOrderStorage(complex<double>*& NextData){

    if(NextData==0){
        // copy into ordered storage
        vector<complex<double> >* stor = new vector<complex<double> >();
        for(Coefficients * l=const_cast<Coefficients*>(firstLeaf());l!=0;l=l->nextLeaf()){
            for(unsigned int k=0;k<l->size();k++)stor->push_back(l->floorData()[k]);
        }
        delete _centralStorage[this];
        _centralStorage[this]=stor;
        NextData=stor->data();
    } else if(_centralStorage.count(this)){
        delete _centralStorage[this];
        _centralStorage.erase(this);
    }

    _cData=NextData;
    if(isLeaf()){
        NextData+=idx()->size();
    } else {
        for(size_t k=0;k<childSize();k++)
            child(k)->treeOrderStorage(NextData);
    }
}

void Coefficients::setOrderedData(std::complex<double> * CData){

    if(_centralStorage.count(this))
    {
        if(CData!=0)DEVABORT("found new storage with CData already set");
        CData=storageData();
        root()->unsetOrderedData(depth()); // cannot be contiguously ordered above
    }
    _cData=CData;

    for(size_t k=0;k<childSize();k++){
        child(k)->setOrderedData(CData);
        if(CData!=0)CData+=child(k)->size(); // advance to beginning of next child
    }
}

void Coefficients::unsetOrderedData(int Level){
    if(isLeaf() or Level<=0)return;
    _cData=0;
    for(size_t k=0;k<childSize();k++)child(k)->unsetOrderedData(Level-1);
}

std::complex<double> *Coefficients::floorData() const {
    if(not isLeaf())return 0;
    return const_cast<Coefficients*>(this)->_cData;
}
void Coefficients::setFloorData(complex<double>* Data) {
    if(not isLeaf())DEVABORT("floor Data only on leaf");
    _cData=Data;
}

std::complex<double> *Coefficients::orderedData() const {
    return const_cast<Coefficients*>(this)->_cData;
}
std::complex<double> *Coefficients::anyData() const {
    if(not notNull())return 0;
    if(_cData!=0)return const_cast<Coefficients*>(this)->_cData;
    auto p=_centralStorage.find(this);
    if(p!=_centralStorage.end())return p->second->data();
    return 0;
}

bool Coefficients::hasFloorData() const {
    if(not isLeaf())return false;
    return _cData!=0;
}

std::complex<double>* Coefficients::data(){
    if(isLeaf() and hasFloorData())return floorData();
    if (storageSize()!=0)return nodeStorage()->data();
    return 0;
}

const std::complex<double>* Coefficients::data() const {
    if(isLeaf() and hasFloorData())return floorData();
    else if (storageSize()!=0)return nodeStorage()->data();
    return 0;
}

TIMERRECURSIVE(copy,)
Coefficients::Coefficients(const Coefficients &Other, std::complex<double>* CDataBegin)
    : Coefficients() {
    setIdx(Other.idx());

    if(Other.storageSize()>0){
        storageAssign(*Other.nodeStorage());
        CDataBegin=nodeStorage()->data();
    }

    if(Other.floorData())_cData=CDataBegin;

    for(unsigned int n = 0; n < Other.childSize(); n++){
        childAdd(new Coefficients(*Other.child(n),CDataBegin));
        if(CDataBegin!=0)CDataBegin+=childBack()->size();
    }
}

void Coefficients::nodeCopy(const Coefficients *Other, bool View)
{
    // CAUTION: does not work correctly when isolated floor
    //          changes were found to be very delicate
    setIdx(Other->idx());
    if(isLeaf())_cData=Other->_cData;
    else        _cData=0;

    if(View){
        makeView();
    } else {
        if(Other->isView()){
            if(idx()->hasFloor())storageAssign(std::vector<std::complex<double>>(Other->_cData,Other->_cData+Other->size()));
        }
        else {
            if(Other->storageSize()>0)storageAssign(*Other->nodeStorage());
            if(Other->hasFloorData()){
                setFloorData(storageData()+(Other->floorData()-Other->storageData()));
            }
        }
    }
}


bool Coefficients::nodeEquivalent(const Coefficients * Other) const {return idx()->nodeEquivalent(Other->idx());}
bool Coefficients::nodeEmpty() const{return not hasFloorData();}

Coefficients::Coefficients(const Index *Idx,const std::vector<std::complex<double>> & Vals):
    Coefficients(Idx){
    if(size()!=Vals.size())
        DEVABORT(Sstr+"initialization size"+Vals.size()+" does not match Coefficients size"+size());
    for(size_t k=0;k<size();k++)data()[k]=Vals[k];
}
Coefficients::Coefficients(const Index *Idx, complex<double> Val)
    :Coefficients(Idx,Val,0){
    treeOrderStorage();
}
Coefficients::Coefficients(const Index *Idx, complex<double> Val,complex<double>*CDataBegin)
    : Coefficients()
{
    parentRef()=0;

    setIdx(Idx);

    if(Idx==0)return; // dummy Coefficients

    // recursively set up coefficient hierarchy
    if (Idx->hasFloor() or Idx->isLeaf()) {
        storageAssign(idx()->size(),Val);
        setOrderedData();
        CDataBegin=floorData();
    }
    else if(Idx->isRoot() or not Idx->parent()->hasFloor()){
        // while above floor, continue descend
        for (unsigned int n=0;n<Idx->childSize();n++){
            childAdd(new Coefficients(Idx->child(n),Val,CDataBegin));
        }
    }
    if(CDataBegin!=0 and idx()->posInFloor()!=Index::npos)
        setFloorData(CDataBegin+idx()->posInFloor());

#ifdef _PARALLEL_
    Parallel::addLeaf(this);
#endif
    if(Idx->parent()==0)treeOrderStorage();
}

Coefficients::Coefficients(std::string File, const Index *&NewIndex /** pointer to "own" the Coefficient's Index */)
    :Coefficients(){
    ifstream istream(File.c_str(),(ios_base::openmode) ios::beg|ios::binary);
    if(not istream.is_open())ABORT("could not find Coefficients file \""+File+"\"");

    NewIndex=new Index(istream,true);
    reset(NewIndex);
    read(istream,false);
}

Coefficients::Coefficients(int FloorDepth, const Index *Idx, complex<double> Val, complex<double>*CDataBegin)
    :Coefficients()
{
    parentRef()=0;
    setIdx(Idx);

    if(Idx==0)return; // dummy Coefficients

    // recursively set up coefficient hierarchy
    if (FloorDepth==0 or (CDataBegin==0 and FloorDepth<0 and Idx->hasFloor())) {
        storageAssign(idx()->size(),Val);
        setOrderedData();
        CDataBegin=floorData();
    }
    else if(Idx->isRoot() or not Idx->parent()->hasFloor()){
        // while above floor, continue descend
        for (unsigned int n=0;n<Idx->childSize();n++){
            childAdd(new Coefficients(FloorDepth-1,Idx->child(n),Val,CDataBegin));
        }
    }
}

//}
void Coefficients::cleanUp(){
    for(auto p: _centralStorage)delete p.second;
    _centralStorage.clear();
}

Coefficients::~Coefficients() {

    // if it was labelled, unregister label (in case pointer is re-assigned to new object)
    auto q=_labels.find(this);
    if(q!=_labels.end())_labels.erase(q);

    auto p=_centralStorage.find(this);
    if(p!=_centralStorage.end()){
        delete p->second;
        _centralStorage.erase(p);
    }

    // remove from list of floor views
    auto f=_viewFloor.find(this);
    if(f!=_viewFloor.end())_viewFloor.erase(f);

    // explicit handling of segfault
    //    signal(SIGSEGV, previous_sigsegv_function);
}

unsigned long Coefficients::size() const {
    return idx()->sizeStored();
}

complex<double> Coefficients::scalarProduct(const Coefficients &RightHandVector) const{
    return idx()->localOverlap()->matrixElement(*this,RightHandVector);
}

complex< double > Coefficients::floorInnerProduct(const Coefficients *ket, bool pseudoScalar) const{
    if(not floorData() or not ket->floorData())return 0.;
    if(pseudoScalar) return (Map<VectorXcd>(floorData(),size())).transpose()*Map<VectorXcd>(ket->floorData(),size());
    else             return (Map<VectorXcd>(floorData(),size())).adjoint()  *Map<VectorXcd>(ket->floorData(),size());
}

complex< double > Coefficients::innerProduct(const Coefficients *ket, bool pseudoScalar) const {
    if (floorData() and ket->floorData()) {
        complex<double> result = floorInnerProduct(ket,pseudoScalar);
        return result;
    }
    else{
        complex< double > result = 0.0;
        if(childSize()!=ket->childSize()) ABORT("Size mismatch");
        for (unsigned int k=0; k<childSize(); k++){
            if(child(k)->notNull() and ket->child(k)->notNull())
                result += child(k)->innerProduct(ket->child(k),pseudoScalar);
        }
        return result;
    }
}

complex< double > Coefficients::innerProductUnscaled(const Coefficients *ket) const {
    if (floorData() and ket->floorData()) {
        // must be floor - need position info
        for(const Index* ix=idx();ix!=0;ix=ix->descend())
            if(idx()->basis()->integrable() and ix->basis()->isAbsorptive())return 0.;

        complex<double> result = floorInnerProduct(ket,false);
        return result;
    }
    else{
        complex< double > result = 0.0;
        if(childSize()!=ket->childSize()) ABORT("Size mismatch: "+tools::str(childSize())+" != "+tools::str(ket->childSize()));
        for (unsigned int k=0; k<childSize(); k++){
            if(child(k)->notNull() and ket->child(k)->notNull())
                result += child(k)->innerProductUnscaled(ket->child(k));
        }
        return result;
    }
}

void Coefficients::makeContinuous(double Scale){

    // NOTE: finite element levels must not be in floor
    if(idx()->hasFloor())return;

    // make sure all lower levels are continuous
    for(unsigned int k=0;k<childSize();k++)child(k)->makeContinuous(Scale);

    if(idx()->continuity()!=Index::npos){

        // average local margins
        for(unsigned int k=1;k<childSize();k++){
            child(k-1)->averageMargin(child(k),idx()->continuity(),Scale);
        }

        // periodic boundaries
        if(child(0)->idx()->basis()->isPeriodic()){
            childBack()->averageMargin(child(0),idx()->continuity(),Scale);
        }
    }
}

void Coefficients::averageMargin(Coefficients* Upper, unsigned int Level, double Scale){

    if(idx()->hasFloor()){
        // apply at floor level
        averageFloorMargin(idx(),floorData(),Upper->idx(),Upper->floorData(),Level,Scale);
    } else {
        if(depth()<Level){
            // descend in Coefficients hierarcy
            for(unsigned int k=0;k<min(childSize(),Upper->childSize());k++)child(k)->averageMargin(Upper->child(k),Level,Scale);
        }
    }
}

void Coefficients::averageFloorMargin(const Index *ILow, std::complex<double> *CLow, const Index *IUpp, std::complex<double> *CUpp, unsigned int Level,const double Scale){
    if(ILow->depth()==Level){
        unsigned int kLow=ILow->basis()->integrable()->upperMargin(),kUpp=IUpp->basis()->integrable()->lowerMargin();
        complex<double>*low;
        complex<double>*upp;
        if(ILow->isBottom()){
            low=CLow+ILow->posInFloor()+kLow;
            upp=CUpp+IUpp->posInFloor()+kUpp;
            *upp=*low=Scale*0.5*(*low+*upp);
        } else {
            low=CLow+ILow->child(kLow)->posInFloor();
            upp=CUpp+IUpp->child(kUpp)->posInFloor();
            for(unsigned int k=0;k<ILow->child(kLow)->size();k++,upp++,low++)*upp=*low=Scale*0.5*(*low+*upp);
        }
    }
    else {
        for(unsigned int k=0;k<ILow->childSize();k++)
            averageFloorMargin(ILow->child(k),CLow,IUpp->child(k),CUpp,Level,Scale);
    }
}

Coefficients* Coefficients::retrieve(const Index *Idx){
    if(idx()==Idx)return const_cast<Coefficients*>(this);
    if(Idx->depth()<idx()->depth())ABORT("below present - cannot retrieve coefficient");
    for(unsigned int k=0;k<childSize();k++){
        if(child(k)->idx()->nSibling()==int(Idx->index()[idx()->depth()]))
            return child(k)->retrieve(Idx);
    }
    ABORT("index not in coefficient: "+tools::str(Idx->index()));
    return 0;
}

void Coefficients::reset(const Index *Idx){
    // delete all branches
    clear();

    // Coefficients with desired structure
    Coefficients rst(Idx);
    //    rst.treeOrderStorage(); // until this is default

    // swap into present
    nodeCopy(&rst,false);
    for(size_t k=0;k<rst.childSize();k++){
        childAdd(rst.child(k));
        rst.child(k)=0; // hide data from being delete'd
    }
    // the solution above probably now fixes the problem correctly
    if(rst.isLeaf()){
        storageAssign(size(),0.);
        setFloorData(storageData());
    } else {
        complex<double>* NextData(0);
        treeOrderStorage(NextData);
    }
}

bool Coefficients::isZero(double eps) const {
    if(notNull()){
        if(anyData())
            return Map<VectorXcd>(anyData(),size()).isZero(eps);
        else
            for (unsigned int k=0; k<childSize(); k++)if(not child(k)->isZero(eps)) return false;
    }
    return true;
}

bool Coefficients::isNan() const {
    if(anyData()){
        for(size_t k=0;k<size();k++)
            if(std::isnan(anyData()[k].real()) or std::isnan(anyData()->imag()))return true;
    }
    else
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull() and child(k)->isNan()) return true;
    return false;
}

void Coefficients::purgeNearZeros(double Eps){
    if(isZero(Eps))setToZero();
    if(floorData()){
        for(unsigned int k=0;k<size();k++)
            if(abs(floorData()[k])<Eps)floorData()[k]=0.;
    }
    else
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull())child(k)->purgeNearZeros(Eps);

}

void Coefficients::setToZero() {
    if (anyData())
        //        memset(anyData(),0,size()*sizeof(*floorData()));
        for(std::complex<double>* a=anyData();a!=anyData()+size();a++)*a=0.;
    else
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull())child(k)->setToZero();
}


void Coefficients::setToConstant(complex<double> Val){
    if (anyData()!=0)
        for(complex<double>*y=anyData();y!=anyData()+size();y++)*y=Val;
    else
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull())child(k)->setToConstant(Val);
}

void Coefficients::setToRandom(){
    if (anyData())
        for(complex<double>*y=anyData();y!=anyData()+size();y++)*y=complex<double>(drand48(),drand48());
    else
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull())child(k)->setToRandom();
}


std::string Coefficients::setToFunctionGrid(std::string Function){
    if(dynamic_cast<CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");
    std::string fail,func;
    std::vector<double> mesh;
    if(idx()->basis()->grid())mesh=idx()->basis()->grid()->mesh();
    if(isLeaf()){
        func=tools::substringReplaceAll(Function,idx()->axisName(),std::string("(Q)"));
        Algebra alg(func);
        if(not alg.isAlgebra())return func;
        for(size_t k=0;k<size();k++)data()[k]=mesh.size()>0?alg.val(mesh[k]):alg.val(0.);
    }
    else {
        for(size_t k=0;k<childSize();k++){
            // replace axisName with mesh value
            func=mesh.size()>0?
                        tools::substringReplaceAll(Function,idx()->axisName(),std::string("("+tools::str(mesh[k])+")"))
                      :Function;
            fail=child(k)->setToFunctionGrid(func);
            if(fail!="")break;
        }
    }
    return fail;
}

/// Coefficients for Function in present basis: Function ~ |i>c[i]
///
/// Function can be any algebraic combination of axis names
/// <br> as e.g. pow[2](Eta)*Rn*exp(-Rn)
/// <br> see class Algebra for available functions
void Coefficients::setToFunction(string Function){
    if(dynamic_cast<CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");

    if(idx()->isHybrid())DEVABORT("must not use with Hybrid basis - set on branch instead, hiearchy: "+idx()->hierarchy());

    PrintOutput::DEVwarning("setToFunction is poorly tested - verify result");

    // get all integrable axes
    std::vector<std::string> axes;
    for(const Coefficients* c=this;c!=0;c=c->nodeNext()){
        if(c->idx()->basis()->integrable() and
                std::find(axes.begin(),axes.end(),c->idx()->axisName())==axes.end()){
            axes.push_back(c->idx()->axisName());
        }
    }
    // get quadratur grid discretiation and vector
    DiscretizationGrid dGrid(idx(),axes);
    Coefficients g(dGrid.idx());

    std::string fail=g.setToFunctionGrid(Function);
    if(fail!="")ABORT("incorrectly expanded function "+Function+" to "+fail+"check algebra, try place brackets around axis names");
    dGrid.mapToParent()->apply(1,g,0,*this);


    /*
    // check FFT for function
    vector<complex<double> > func;
    double L=40.;
    int ptsIn=40;
    for(size_t k=0;k<ptsIn;k++)func.push_back(alg.val(k*L/ptsIn));
    Fft fft(func.size(),true);

    vector<complex<double> > ft=func;
    fft.transform(ft);

    Fft bft(func.size(),false);

    ofstream out,ftrans;
    ftrans.open("trans");
    ftrans<<"# k ft"<<endl;
    for(size_t k=0;k<ptsIn;k++)ftrans<<k<<", "<<std::norm(ft[k])<<endl;

    vector<complex<double> > bt=ft;
    bft.transform(bt);

    out.open("fourier");
    out<<"# axis plot "+Function+" func bt ft"<<endl;

    int pts=ptsIn*1.7;
    for(size_t l=0;l<pts;l++){
        double x=l*L/pts;
        complex<double> p(0.);
        for(size_t k=0;k<ptsIn/2;k++)p+=exp(complex<double>(0.,2*math::pi/L*k*x))*ft[k];
        for(size_t k=ptsIn/2;k<ptsIn;k++)p+=exp(complex<double>(0.,2*math::pi/L*(k-ptsIn)*x))*ft[k];
        out<<x<<", "<<std::norm(p)<<", "<<std::norm(alg.val(x))<<endl;
    }
*/
}

void Coefficients::conjugate(){
    if (anyData())
        for(complex<double>*y=anyData();y!=anyData()+size();y++)*y=std::conj(*y);
    else
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull())child(k)->conjugate();
}

Coefficients *Coefficients::firstFloor() const {
    if(dynamic_cast<const CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");
    if(hasFloorData() or childSize()==0)return const_cast<Coefficients*>(this);
    return child(0)->firstFloor();
}

void Coefficients::pointerToC(vector<complex<double> *> & pC) {
    if(dynamic_cast<CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");

    // top level:
    if(idx()->parent()==0 and pC.size()>0){
        // check whether pC contains right indices
        const Coefficients * leaf=firstFloor();
        while(not leaf->hasFloorData()){
            leaf=leaf->nextLeaf();
            if(leaf==0)ABORT("coefficients do not contain any floors");
        }
        if(pC[0]==leaf->floorData())return; //
        pC.clear();
    }

    if (hasFloorData())
        for (unsigned int n=0;n<size();n++)
            pC.push_back(&(floorData()[n]));

    for (unsigned int n=0;n<childSize();n++)
        child(n)->pointerToC(pC);
}

Coefficients& Coefficients::operator=(const Coefficients& rhs)
{
    if (this == &rhs) return *this;
    if(idx()==0){
        reset(rhs.idx());
    }
    if((orderedData()!=0) and (rhs.orderedData()!=0)){
        if(size()!=rhs.size())
            DEVABORT(Sstr+"vector sizes differ - cannot assign:"+size()+"<--"+rhs.size()+"\nrhs\n"+rhs.idx()->str()+"\nlhs\n"+idx()->str());
        if(size()==0)PrintOutput::DEVwarning("empty assign");
        if(orderedData()!=rhs.orderedData())memcpy(orderedData(),rhs.orderedData(),size()*sizeof(*orderedData()));
    }
    else
        for (unsigned int k=0; k<rhs.childSize(); k++)
            if(child(k)->notNull())*child(k)=*rhs.child(k);
    return *this;
}

Coefficients &Coefficients::cwiseMultiply(const Coefficients& B)
{
    if((orderedData()!=0) and (B.orderedData()!=0)){
        for(size_t k=0;k<size();k++)orderedData()[k]*=B.orderedData()[k];
    }
    else
        for (unsigned int k=0; k<B.childSize(); k++)
            if(child(k)->notNull() and B.child(k)->notNull())child(k)->cwiseMultiply(*B.child(k));
    return *this;
}

Coefficients & Coefficients::cwiseDivide(const Coefficients& B)
{
    //    if((floorData()!=0) and (B.floorData()!=0))
    if((orderedData()!=0) and (B.orderedData()!=0))
        for(size_t k=0;k<size();k++)orderedData()[k]/=B.orderedData()[k];
    else
        for (unsigned int k=0; k<B.childSize(); k++)
            if(child(k)->notNull() and B.child(k)->notNull())child(k)->cwiseDivide(*B.child(k));
    return *this;
}

Coefficients & Coefficients::cwiseRelativeError(const Coefficients& B, double Bound)
{
    if((orderedData()!=0) and (B.orderedData()!=0))
        for(size_t k=0;k<size();k++)orderedData()[k]=
                2.*(orderedData()[k]-B.orderedData()[k])/std::max(std::abs(orderedData()[k]+B.orderedData()[k]),Bound);
    else
        for (unsigned int k=0; k<B.childSize(); k++)
            if(child(k)->notNull() and B.child(k)->notNull())child(k)->cwiseRelativeError(*B.child(k),Bound);
    return *this;
}

Coefficients& Coefficients::operator+=(const Coefficients& rhs)
{
    if (this == &rhs)ABORT("no aliasing in += allowed");
    if((orderedData()!=0) and (rhs.orderedData()!=0))
        for(complex<double>*y=orderedData(),*x=rhs.orderedData();y!=orderedData()+size();y++,x++)*y+=*x;
    else{
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull() and rhs.child(k)->notNull())*child(k) += *rhs.child(k);
    }
    return *this;
}

Coefficients& Coefficients::operator-=(const Coefficients& rhs)
{
    if (this == &rhs)ABORT("no aliasing in -= allowed");
    if((orderedData()!=0) and (rhs.orderedData()!=0)){
        for(complex<double>*y=orderedData(),*x=rhs.orderedData();y!=orderedData()+size();y++,x++)*y-=*x;
    }
    else{
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull() and rhs.child(k))*child(k)-=*rhs.child(k);
    }
    return *this;
}


Coefficients& Coefficients::operator*=(std::complex< double > A)
{
    if(anyData()!=0)
        for(complex<double>*y=anyData();y!=anyData()+size();y++)*y*=A;
    else
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull())*child(k)*=A;

    return *this;
}

Coefficients & Coefficients::axpy(std::complex<double> A, const Coefficients &X, std::complex<double> B)
{
    if(A==0.){
        if(B==0.)setToZero();
        else if(B!=1.)scale(B);
        return *this;
    }

    if(B==1.)
        axpy(A,&X);
    else {
        if ((orderedData()!=0) and (X.orderedData()!=0)){
            if(A==1.)for(complex<double>*y=orderedData(),*x=X.orderedData();y!=orderedData()+size();y++,x++)(*y*=B)+=(*x);
            else     for(complex<double>*y=orderedData(),*x=X.orderedData();y!=orderedData()+size();y++,x++)(*y*=B)+=(*x)*A;
        }
        else
            for (unsigned int k=0; k<childSize(); k++)
                if(child(k)->notNull() and X.child(k)->notNull())child(k)->axpy(A,*X.child(k),B);
    }
    return *this;
}

void Coefficients::axpy(std::complex< double > A, const Coefficients * X)
{
    if(A==0.)return;
    if ((orderedData()!=0) and (X->orderedData()!=0))
        for(complex<double>*y=orderedData(),*x=X->orderedData();y!=orderedData()+size();y++,x++)*y+=(*x)*A;
    else
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull() and X->child(k)->notNull())child(k)->axpy(A,X->child(k));
}

void Coefficients::axpy(const Eigen::MatrixXcd &Amat, const Coefficients *X){
    if(Amat.rows()!=childSize() or Amat.cols()!=X->childSize())
        ABORT(Str("matrix does not match coefs")
              +Amat.rows()+"X"+Amat.cols()+"vs."+childSize()+"X"+X->childSize());
    for(size_t i=0;i<childSize();i++)
        for(size_t j=0;j<X->childSize();j++)
            if(Amat(i,j)!=0.)
                if(child(i)->notNull() and X->child(j)->notNull())child(i)->axpy(Amat(i,j),X->child(j));
}

void Coefficients::scale(const std::complex< double > A)
{
    if(A==1.)return;
    if(A==0.){
        setToZero();
        return;
    }
    if (anyData())
        for(complex<double>*y=anyData();y!=anyData()+size();y++)*y*=A;
    else
        for (unsigned int k=0; k<childSize(); k++)
            if(child(k)->notNull())child(k)->scale(A);
}

std::complex<double> Coefficients::cMaxNorm() const {
    std::complex<double> cMax(0.);
    if(anyData()){
        int kPos=std::max_element(anyData(),anyData()+size(),
                                  [](std::complex<double>a,std::complex<double>b){return std::norm(a)<std::norm(b);})-anyData();
        cMax=anyData()[kPos];
    }
    else
        for(size_t k=0;k<childSize();k++){
            if(child(k)->notNull()){
                std::complex<double> cChild=child(k)->cMaxNorm();
                if(std::abs(cMax)<std::abs(cChild))cMax=cChild;
            }
        }
    return cMax;
}


std::vector<std::complex<double> > *Coefficients::nodeStorage() const {
    unordered_map<const Coefficients*,vector<complex<double> >* >::iterator p=_centralStorage.find(this);
    if(p==_centralStorage.end())return 0;
    return p->second;
}

void Coefficients::storageAssign(int Size, std::complex<double> Val) const {
    storageAssign(std::vector<std::complex<double> >(Size,Val));
}
void Coefficients::storageAssign(const std::vector<std::complex<double> > & Stor) const {
    unordered_map<const Coefficients*,vector<complex<double> >* >::iterator p=_centralStorage.find(this);
    if(p!=_centralStorage.end()){
        delete p->second;
        p->second=new vector<complex<double> >(Stor);
        // previous pointer looses meaning, new storage is not guaranteed to be ordered
        const_cast<Coefficients*>(this)->_cData=0;
    }
    else
        _centralStorage[this]=new vector<complex<double> >(Stor);

    if(Stor.size()==1)
        const_cast<Coefficients*>(this)->_cData=_centralStorage[this]->data();
}

void Coefficients::nodeStorageClear() const {
    unordered_map<const Coefficients*,vector<complex<double> >* >::const_iterator p=_centralStorage.find(this);
    if(p!=_centralStorage.end()){
        p->second->clear();
        delete p->second;
        _centralStorage.erase(p);
    }
}
unsigned int Coefficients::storageSize() const {
    unordered_map<const Coefficients*,vector<complex<double> >* >::const_iterator p=_centralStorage.find(this);
    if(p==_centralStorage.end())return 0;
    return p->second->size();
}

double Coefficients::norm() const {
    double nrm=-1.;
    switch (normType){
    case 0:
        if (orderedData())
            for(size_t k=0;k<size();k++)nrm=max(nrm,max(abs(orderedData()[k].real()),abs(anyData()[k].imag())));
        else
            for(unsigned int k=0;k<childSize();k++)
                if(child(k)->notNull())nrm=max(nrm,child(k)->norm());
        break;
    case 2:
        nrm=0.;
        if(orderedData())
            for(size_t k=0;k<size();k++)nrm+=std::norm(orderedData()[k]);
        else
            for(unsigned int k=0;k<childSize();k++)
                if(child(k)->notNull())nrm+=child(k)->norm();
        break;
    default:
        ABORT("only norm types 0(=infty norm) and 2(=square[L2-norm]) defined, is: "+tools::str(normType));
    }
    return nrm;
}

double Coefficients::pNorm(int P) const {
    double nrm=-1.;
    switch (P){
    case 0:
        if (orderedData())
            for(size_t k=0;k<size();k++)nrm=max(nrm,std::abs(orderedData()[k]));
        else
            for(unsigned int k=0;k<childSize();k++)
                if(child(k)->notNull())nrm=max(nrm,child(k)->pNorm(P));
        break;
    case 1:
        nrm=0.;
        if(orderedData())
            for(size_t k=0;k<size();k++)nrm+=std::abs(orderedData()[k]);
        else
            for(unsigned int k=0;k<childSize();k++)
                if(child(k)->notNull())nrm+=child(k)->pNorm(P);
        break;
    case 2:
        nrm=0.;
        if(orderedData())
            for(size_t k=0;k<size();k++)nrm+=std::norm(orderedData()[k]);
        else
            for(unsigned int k=0;k<childSize();k++)
                if(child(k)->notNull())nrm+=child(k)->pNorm(P);
        break;
    default:
        ABORT("only P=0(=infty norm), 1(=1-norm) and 2(=square[L2-norm]) defined, is: "+tools::str(P));
    }
    return nrm;
}


std::vector<std::complex<double>> Coefficients::values() const {
    if(dynamic_cast<const CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");
    std::vector<std::complex<double>> v;
    if(orderedData())
        for(size_t k=0;k<size();k++)v.push_back(orderedData()[k]);
    else{
        for(size_t k=0;k<childSize();k++){
            std::vector<std::complex<double>> cv(child(k)->values());
            v.insert(v.end(),cv.begin(),cv.end());
        }
    }
    return v;
}

// this should go into an auxiliary namespace/class
void Coefficients::plot(std::string PltFile){
    if(dynamic_cast<const CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");
    if(size()==0)return;
    Plot plt(idx(),ReadInput::main);
    if(plt.isEmpty())return;
    if(PltFile==""){
        size_t cnt=0;
        PltFile=ReadInput::main.outputTopDir()+"C_00";
        while(folder::exists(PltFile) and cnt<99)PltFile=PltFile.substr(0,PltFile.length()-2)+tools::str(cnt,2,'0');
    }
    plt.plot(*this,PltFile);
    PrintOutput::message("Coeffients plot on "+PltFile);
}
void Coefficients::makeView() const{
    if(idx()->hasFloor())
        _viewFloor.insert(this);
}

/// Precision=0..character symbols,
/// =floorPointers..pointers,
/// =n..digits,=100+n...n digits of abs(C),
/// =200+n...n digits of arg(C)/pi */
string Coefficients::strNode(int Precision) const{
    if(dynamic_cast<const CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");

    Str s("","");
    if(idx()==0)ABORT("idx() not assigned");
    s+=" "+idx()->strNode();
    if(childSize()>0 or hasFloorData()){
        if(          orderedData()!=0)      s=s+" pData: "+orderedData();
        else if(_centralStorage.count(this))s=s+" pStor: "+_centralStorage[this]->data();
        else                                s=s+" noData";
        if(isView())s=s+":V";
        s=s+" nrm="+tools::str(Coefficients::norm(),2)+" siz="+size();
        if(childSize()>0)
            s=s+" sub="+childSize();
        else {
            if(MPIwrapper::Size()>1){
                s=s+" <"+Parallel::owner(idx());
                if(dynamic_cast<const CoefficientsLocal*>(root()))s=s+"@"+idx()->posIndex(root()->idx());
                s=s+">";
            }
            s+="\t";
            if(Precision!=Tree_ptrsSizes){
                for(unsigned int n=0;n<size();n++){
                    if(floorData()!=0){
                        if (Precision==Tree_defaultKind)s=s+tools::zero(floorData()[n]);
                        else if(Precision==floorPointers)s=s+(floorData()+n)+" ";
                        else if (Precision==0)s=s+tools::zero(floorData()[n]);
                        else if(Precision<100)s=s+tools::str(floorData()[n],Precision)+" ";
                        else if(Precision<200)s=s+tools::str(abs(floorData()[n]),Precision-100)+" ";
                        else if(Precision<300)s=s+tools::str(std::arg(floorData()[n])/math::pi,Precision-200)+" ";
                        else if(Precision<400)s=s+tools::str(std::complex<double>(std::abs(floorData()[n]),std::arg(floorData()[n])/math::pi),Precision-300)+" ";
                    } else {
                        s=s+"(no data)";
                    }

                }
            }
        }
    }
    else {
        s+=" (empty leaf)";
    }
    return std::move(s);
}

/// recursively write to binary: ofstream stream(filename,.c_str(),(ios_base::openmode) ios::beg|ios::binary)
void Coefficients::write(ofstream &Stream, bool Header, string Kind) const {
    if(dynamic_cast<const CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");

    if(idx()->parent()==0){
        //        if(idx()->str().find("subset")!=string::npos)ABORT(idx()->str()+"\nmust not write threaded");
        if(not Stream.is_open())ABORT("need to open output stream before writing coefficients");
        if(Header){
            Stream.seekp(ios_base::beg); // header is at beginning of file
            if(Kind=="IndexStructure")idx()->writeStructure(Stream);
            else if(Kind=="IndexFull")idx()->write(Stream);
            else ABORT("Undefined Kind="+Kind+", admissibel IndexStructur|IndexFull");
        }
    }
    if(isLeaf() and hasFloorData())tools::write(Stream,_cData,size());
    for (unsigned int n=0;n<childSize();n++)child(n)->write(Stream,false);

    if(idx()->parent()==0){
        Stream.flush(); // when exiting from top level, flush data
    }
}

/// recursively write to ascii ofstream stream(filename,.c_str(),(ios_base::openmode) ios::beg|ios::binary)
void Coefficients::print(ofstream &stream) const {
    if(dynamic_cast<const CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");
    if(idx()->parent()==0 and not stream.is_open())ABORT("need to open output stream before writing coefficients");
    if(isLeaf())
        for(unsigned int k=0;k<size();k++)stream<<", "<<abs(*(_cData+k))<<", "<<arg(*(_cData+k));
    for (unsigned int n=0;n<childSize();n++)child(n)->print(stream);

    if(idx()->parent()==0){
        stream<<endl;
    }

}

/// recursively read from binary file
bool Coefficients::read(ifstream &stream, bool header){
    if(dynamic_cast<const CoefficientsSparse*>(this))DEVABORT("cannot use CoefficientsSparse here");
    if(header){ // check header information
        stream.seekg(ios_base::beg); // header is at beginning of file
        if(not idx()->compatibleFile(stream))return false;
    }
    if(isLeaf()){
        tools::read(stream,_cData,size());
    }
    for (unsigned int n=0;n<childSize();n++)
        if(not child(n)->read(stream,false))return false;
    return stream.good();

}

void Coefficients::print(const std::vector<Coefficients *> Coeff, string Text){
    if(Coeff.size()==0){
        PrintOutput::warning("Coefficients "+Text+": emtpy");
        return;
    }
    UseMatrix vec(Coeff[0]->size(),Coeff.size());
    vector<complex<double>* > pC(Coeff[0]->size());
    for(unsigned int k=0;k<Coeff.size();k++){
        if(dynamic_cast<const CoefficientsSparse*>(Coeff[k]))DEVABORT("cannot use CoefficientsSparse here");
        Coeff[k]->pointerToC(pC);
        for(unsigned int l=0;l<pC.size();l++)vec(l,k)=*pC[l];
    }
    vec.print(Text);
}


void Coefficients::examplePermute(){

    //------------------------------------------
    // interchange the indices
    //----------------------------------------

    // create index view with floor moved to bottom
    vector<unsigned int> perm;
    for(unsigned int k=0;k<idx()->firstLeaf()->depth();k++)perm.push_back(k);
    std::swap(perm[perm.size()-2],perm[perm.size()-1]);
    cout<<"permutation "<<tools::str(perm)<<endl;

    cout<<"source "<<idx()->str()<<endl;

    Index iSrc0(*idx());
    iSrc0.resetFloor(iSrc0.firstLeaf()->depth());
    cout<<"floor lowered\n"<<iSrc0.str()<<endl;
    Index pSrc0;
    iSrc0.permute(perm,pSrc0);
    cout<<"permuted\n"<<pSrc0.str()<<endl;

    // create index with lower levels merged into floor
    Index iTarget(pSrc0);
    iTarget.resetFloor(iTarget.firstLeaf()->depth()-2);
    cout<<"target index\n"<<iTarget.str()<<endl;


    //----------------------------------------
    // create coefficients and permuted view
    //----------------------------------------

    // create source coefficients
    Coefficients cSource(idx());

    // create target coefficients
    Coefficients cTarget(&iTarget);

    // create views on source and target with floor lowered
    Coefficients vSource(&iSrc0,&cSource);
    Coefficients vTarget(&pSrc0,&cTarget);

    // create view on permuted source coefficients
    Coefficients pSource;
    vSource.permute(perm,pSource);
    Coefficients wTarg0(&pSrc0,&pSource);

    //-------------------------------------------------------
    // assign values to source (and to target, for checking)
    //-------------------------------------------------------
    vector<complex<double>*> pC;
    cSource.pointerToC(pC);
    for(unsigned int k=0;k<pC.size();k++)*(pC[k])=k;

    cTarget.pointerToC(pC);
    for(unsigned int k=0;k<pC.size();k++)*(pC[k])=k*0.1;

    cout<<"source\n"<<cSource.str(2)<<endl;
    cout<<"view on permuted sources coefficients \n"<<wTarg0.str(2)<<endl;

    // copy view on permuted source into target view
    vTarget=wTarg0;

    cout<<"viewTarget \n"<<vTarget.str(2)<<endl;
    cout<<"target \n"<<cTarget.str(2)<<endl;

    PrintOutput::title("permutation demo done");
    exit(0);
}

Coefficients::Coefficients(const Index *I, Coefficients* Origin):Coefficients(){
    if(not Origin)DEVABORT("Coefficients(I,0) called - for initializing to 0. call Coefficients(I,0.)");
    setIdx(I);
    if(Origin->hasFloorData()){
        //CHECK _cData / hasFloorData() usage
        if(idx()->hasFloor()){
            _cData=Origin->_cData;
        }
        else if(Origin->idx()->isLeaf()){
            if(not Origin->idx()->isBottom())DEVABORT("this should not happen: origin index ends, but is not bottom?");
            for(size_t k=0;k<idx()->childSize();k++){
                // construct leaf view - should be cast into function or constructor
                childAdd(new Coefficients(idx()->child(k)));
                childBack()->setFloorData(Origin->orderedData()+k);
                childBack()->makeView();
                _centralStorage.erase(childBack()); // original Coefficients constructor assignes storage - remove
            }
        }
        else {
            for(unsigned int k=0;k<idx()->childSize();k++){
                childAdd(new Coefficients(idx()->child(k),Origin->idx()->child(k),idx()->heightAboveFloor()-1,Origin->floorData()));
            }
        }
    }
    else {
        if(I->hasFloor())ABORT("cannot view from above Coefficients floor");
        if(idx()->childSize()!=Origin->childSize())
            ABORT("Origin Coefficients not compatible with Index: childSize differs");
        for(unsigned int k=0;k<I->childSize();k++)
            childAdd(new Coefficients(I->child(k),Origin->child(k)));
    }
}

Coefficients::Coefficients(Coefficients &Origin, std::function<bool(Coefficients*)> SelectView):Coefficients(){
    setIdx(Origin.idx());
    if(SelectView(&Origin)){
        _cData=Origin._cData;
        makeView();
    }
    else {
        for(unsigned int k=0;k<Origin.childSize();k++)
            if(SelectView(Origin.child(k))){
                childView(Origin.child(k));
            } else {
                childAdd(new Coefficients(*Origin.child(k),SelectView));
            }
    }
}

Coefficients::Coefficients(const Index *I, const Index* IOrigin, unsigned int FloorDepth, std::complex<double> *CDataBegin)
    :Coefficients()
{
    parentRef()=0;
    setIdx(I);
    if(FloorDepth==0){
        unsigned int pos=IOrigin->posInFloor();
        if(pos==Index::npos)ABORT("IOrigin not in floor"+tools::str(IOrigin->index())+"\n"+IOrigin->str()+"\nRoot\n"+IOrigin->root()->str());
        _cData=CDataBegin+IOrigin->posInFloor();
    }
    else {
        for(unsigned int k=0;k<I->childSize();k++)
            childAdd(new Coefficients(I->child(k),IOrigin->child(k),FloorDepth-1,CDataBegin));
    }
}
