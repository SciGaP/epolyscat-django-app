// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "operatorMap.h"

#include "mpiWrapper.h"
#include "typeinfo"
#include "index.h"
#include "str.h"
#include "basisIntegrable.h"
#include "basisSub.h"
#include "coefficients.h"
#include "coefficientsFloor.h"
#include "coefficientsViewDeep.h"
#include "algebra.h"
#include "parallelOperator.h"
#include "parameters.h"
#include "parallel.h"

#include "basisGridQuad.h"
#include "basisSub.h"
#include "eigenTools.h"

using namespace std;
OperatorMap::~OperatorMap(){
    delete _deepX;delete _deepY;delete _highX;delete _highY;delete _viewX;delete _viewY;delete _iX;delete _iY;delete _tmp;
}

OperatorMap::Storage OperatorMap::Storage::main;

const Eigen::MatrixXcd* OperatorMap::Storage::get(Eigen::MatrixXcd&& from){
    for(auto& m: _matrices){
        if(m->rows() == from.rows() and m->cols() == from.cols()){
            if(m->isApprox(from, 1.e-12)) return m.get();
        }
    }

    _matrices.push_back(std::unique_ptr<Eigen::MatrixXcd>(new Eigen::MatrixXcd(from)));
    return _matrices.back().get();
}

// we should write a BasisMap1d base class with the repective instances
static Eigen::MatrixXcd kFromSurf(const Index* IIndex,const Index* JIndex){
    Eigen::MatrixXcd _map;
    string toAx=IIndex->axisName(),fromAx=JIndex->axisName();
    if(not(toAx.find("k")==0 and fromAx.find("ValDer")==0 and toAx.substr(1)==fromAx.substr(6)))return _map;

    if(toAx.find("Rn")==1){
        string axCnt=toAx.substr(3);
        //<spherBessel[0,rightEta,rightPhi,rightsurfRn](Q)>  + <spherBessel[1,rightEta,rightPhi,rightsurfRn](Q)>
        _map.resize(IIndex->basis()->size(),JIndex->basis()->size());
        // determine surface radius etc.
        double radius=JIndex->basis()->grid()->mesh()[0],lAngle=DBL_MAX;
        for(const Index* idx=IIndex;idx->parent()!=0;idx=idx->parent()){
            if(idx->parent()->axisName()=="Eta"+axCnt){
                lAngle=idx->parent()->basis()->physical(idx->nSibling());
            }
        }
        if(lAngle==DBL_MAX)ABORT("did not find Eta"+axCnt+" in hierarchy"+IIndex->root()->hierarchy());

        string lmRn=Str("","")+","+int(lAngle)+","+radius+"]*"+radius;
        Algebra b0("spherBessel[0"+lmRn);
        Algebra b1("spherBessel[1"+lmRn);
        for(size_t k=0;k<IIndex->basis()->size();k++){
            _map(k,0)=b0.val(IIndex->basis()->grid()->mesh()[k]);
            _map(k,1)=b1.val(IIndex->basis()->grid()->mesh()[k]);
        }

    }
    else if(toAx.find("X")==1 or toAx.find("Y")==1 or toAx.find("Z")==1 ){
        //"<Id><ExpI[0,rightsurfX](Q)>+<Id><ExpI[1,rightsurfX](Q)>"
        _map.resize(IIndex->basis()->size(),JIndex->basis()->size());
        double x=JIndex->basis()->grid()->mesh()[0];
        Algebra b0(Str("ExpI[0,","")+x+"](Q)");
        Algebra b1(Str("ExpI[1,","")+x+"](Q)");
        vector<double> g(BasisGrid::factory(IIndex->basis())->mesh());
        for(size_t k=0;k<g.size();k++){
            _map(k,0)=b0.val(-g[k]);
            _map(k,1)=b1.val(-g[k]);
        }
    }
    _map.purge();
    return _map;
}


bool OperatorMap::isIdentity(double Eps, bool Stochastic) const {
    if(_tmp!=0 and not (_tensor!=0 and _tensor->isIdentity(Eps)))return false;
    for(size_t k=0;k<childSize();k++)
        if(not child(k)->isIdentity(Eps))return false;
    if(not iIndex->treeEquivalent(jIndex))return false;
    if(not Stochastic)return true; // shut up unused variable warning
    return true;
}

OperatorMap::OperatorMap(const Index *IIndex, const Index *JIndex, string Derivative, const std::complex<double> Multiplier, bool EntryLevel)
    :OperatorAbstract("Map",IIndex,JIndex),
      _deepX(0),_deepY(0),_highX(0),_highY(0),_viewX(0),_viewY(0),_iX(0),_iY(0),_tmp(0),_tensor(0)
{
    Eigen::MatrixXcd tensor;

    _iVec=IIndex->nSibling();
    _jVec=JIndex->nSibling();

    // check for derivative wrt to current coordinate axis
    bool jDer=JIndex->axisName()==Derivative and JIndex->continuity()==Index::npos; // derivative, if not continuity level

    // special cases
    if((IIndex->axisName().substr(0,6)=="ValDer") != (JIndex->axisName().substr(0,6)=="ValDer"))
        tensor=kFromSurf(IIndex,JIndex)*Multiplier; // map value and derivative into k-grid
    else
        // standard case
        tensor=basisMap(IIndex->basis(),JIndex->basis(),jDer)*Multiplier;

    if(tensor.size()==0){
        ABORT(Str("cannot map\nroots:"," ")
              +IIndex->root()->hierarchy()+"<----"+JIndex->root()->hierarchy()+"\n"
              +IIndex->index()+IIndex->basis()->str()+"<-//-"+JIndex->index()+JIndex->basis()->str());
    }
    if(iIndex->isBottom()){
        if(tensor.size() != 0) _tensor=Storage::main.get(std::move(tensor));
        return;
    }

    bool iEquiv=iIndex->subEquivalent(),jEquiv=jIndex->subEquivalent();
    Coefficients* _aux=0; // preliminary vector
    if(tensor.nonZeros()>min(tensor.rows(),tensor.cols())){
        if(iEquiv and jEquiv){
            // hard to figure out, but assume lower level is much heavier than tensor multipliciation
            // then we should pick the smaller number of lower level applications
            if(iIndex->childSize()>jIndex->childSize())_aux=new Coefficients(iIndex->heightAboveBottom()-1,iIndex->child(0));
            else                                       _aux=new Coefficients(jIndex->heightAboveBottom()-1,jIndex->child(0));
        }
        else if(iEquiv)_aux=new Coefficients(iIndex->heightAboveBottom()-1,iIndex->child(0));
        else if(jEquiv)_aux=new Coefficients(jIndex->heightAboveBottom()-1,jIndex->child(0));
    }

    if(_aux==0){
        for(int j=0;j<tensor.cols();j++){
            for(int i=0;i<tensor.rows();i++){
                if(tensor(i,j)!=0.){
                    childAdd(new OperatorMap(iIndex->child(i),jIndex->child(j),Derivative,tensor(i,j),false));
                }
            }
        }
        tensor.resize(0, 0);
    }
    else if(_aux->idx()==iIndex->child(0)){
        for(size_t j=0;j<jIndex->childSize();j++){
            childAdd(new OperatorMap(_aux->idx(),jIndex->child(j),Derivative,1.,false));
        }
    }
    else if(_aux->idx()==jIndex->child(0)){
        for(size_t i=0;i<iIndex->childSize();i++){
            childAdd(new OperatorMap(iIndex->child(i),_aux->idx(),Derivative,1.,false));
        }
    }
    else
        DEVABORT("algorithm flawed");


    // terminate OperatorMap, if all lower maps are identity
    size_t k;
    for(k=0;k<childSize();k++)
        if(not child(k)->isIdentity(1.e-12))break;
    if(k==childSize())
        while(0<k--)childErase(k);

    // get _tmp: height matches height of OperatorMap children
    if(_aux!=0){
        if(height()==0)
            delete _aux;
        else if(_aux->height()==height()-1){
            _tmp=_aux;
        }
        else {
            _tmp=new Coefficients(height()-1,_aux->idx());
            delete _aux;
        }
    }
    if(EntryLevel){
        if(iIndex->heightAboveBottom()!=jIndex->heightAboveBottom())DEVABORT(Sstr+"left and right index heights differ, hierarchies:"
                                                                             +iIndex->hierarchy()+jIndex->hierarchy());
        size_t firstFloorDepth = max(iIndex->firstFloor()->depth(),jIndex->firstFloor()->depth());
        if(height()>=firstFloorDepth)
        { // views into the floor
            _deepX=new CoefficientsViewDeep(jIndex,height(),true);
            _deepY=new CoefficientsViewDeep(iIndex,height(),true);
        }
        // contiguous storage above floor
        if(height()<firstFloorDepth)
        {
            _iX=new Index(*jIndex);
            _iY=new Index(*iIndex);
            _iX->resetFloor(height());
            _iY->resetFloor(height()); // get indices with floors raised
            _highX=new Coefficients(_iX); // create matching coefficients
            _highY=new Coefficients(_iY,0.);
            _viewX=new Coefficients(jIndex,_highX);// get view matching the input vector
            _viewY=new Coefficients(iIndex,_highY);
        }
    }

    // Actually store tensor
    if(tensor.size() != 0)_tensor=Storage::main.get(std::move(tensor));
}


Eigen::MatrixXcd OperatorMap::basisMap(const BasisAbstract *IBas, const BasisAbstract *JBas, int JDerivative){

    if(JDerivative and dynamic_cast<const BasisIntegrable*>(JBas)==0)
        ABORT("derivative only for differentiable rhs functions");

    if(IBas==JBas)return Eigen::MatrixXcd::Identity(IBas->size(),JBas->size());

    Eigen::MatrixXcd res=BasisSub::map(IBas,JBas);
    if(res.size()!=0)return res;

    if(*IBas==*JBas)return Eigen::MatrixXcd::Identity(IBas->size(),JBas->size());

    const BasisAbstract *iSup=BasisSub::superBas(IBas), *jSup=BasisSub::superBas(JBas);
    if(*iSup==*jSup){
        std::vector<int> iSubs(BasisSub::subset(IBas));
        std::vector<int> jSubs(BasisSub::subset(JBas));
        Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(iSubs.size(),jSubs.size());
        for(size_t j=0;j<jSubs.size();j++)
            for(size_t i=0;i<iSubs.size();i++)
                if(iSubs[i]==jSubs[j])mat(i,j)=1.;
        return mat;
    }

    // map between grids
    if(IBas->isGrid() and JBas->isGrid()){
        const BasisGrid* gi=BasisGrid::factory(IBas);
        const BasisGrid* gj=BasisGrid::factory(JBas);
        if(gi->mesh()==gj->mesh())return Eigen::MatrixXcd::Identity(IBas->size(),JBas->size());
        else                      return gj->mapInterpolate(gi);
    }

    // maps to or from grid
    UseMatrix mat,dum;
    const BasisGrid* iBas=dynamic_cast<const BasisGrid*>(IBas);
    const BasisGridQuad* jBas=dynamic_cast<const BasisGridQuad*>(JBas);
    if(iBas and not jBas){
        if(JDerivative>1)ABORT("for now, only first derivative");
        UseMatrix points(iBas->size(),1);
        for(size_t k=0;k<iBas->size();k++)points(k)=iBas->mesh()[k];
        if(JBas->integrable())JBas->integrable()->valDer(points,mat,dum,JDerivative);
        else if(JBas->sub())         JBas->sub()->valDer(points,mat,dum,JDerivative);
        else DEVABORT(JBas->str()+" is neither integrable nor subset");
    }
    else if (not iBas and jBas){
        UseMatrix tmat;
        UseMatrix points(jBas->size(),1);
        for(size_t k=0;k<JBas->size();k++)points(k)=jBas->mesh()[k];
        if(IBas->integrable())IBas->integrable()->valDer(points,tmat,dum,JDerivative);
        else if(IBas->sub())         IBas->sub()->valDer(points,tmat,dum,JDerivative);
        else DEVABORT(IBas->str()+" is neither integrable nor subset");
        if(IBas->isAbsorptive())mat=tmat.transpose();
        else                    mat=tmat.adjoint();
        for (size_t k=0;k<jBas->weights().size();k++)mat.col(k)*=jBas->weights()[k];
    }
    return Eigen::MatrixXcd::Map(mat.data(),mat.rows(),mat.cols());
}

void OperatorMap::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    if(Vec.idx()!=jIndex)DEVABORT("rhs index does not match\nVector\n"+Vec.idx()->str()+"\nOperator\n"+jIndex->str());
    if(  Y.idx()!=iIndex)DEVABORT("lhs index does not match");

    Y.scale(B);
    if(A==0.)return;

    if(_deepX!=0){
        Coefficients *Xtem = _deepX->view(const_cast<Coefficients*>(&Vec));
        Coefficients * Ytem = _deepY->view(&Y);
        axpy3(A,Xtem,Ytem);
    }
    else if(_highX!=0){
        *_viewX=Vec;
        _highY->setToZero();
        axpy3(A,_highX,_highY);
        Y+=*_viewY;
    }
    else{
        axpy3(A,&Vec,&Y);
    }
}

void OperatorMap::axpy3(std::complex<double> A, const Coefficients *X, Coefficients *Y) const{

    if(isLeaf()){
        if(_tensor==0){
            Y->axpy(A,X);
        }else{
            Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::Map
                    (Y->data(),_tensor->rows(),iIndex->sizeStored()/_tensor->rows())+=
                    *_tensor*
                    (Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::Map
                     (X->data(),_tensor->cols(),jIndex->sizeStored()/_tensor->cols()))*A;
        }
    }
    else {
        // general block structure
        if(_tmp==0){
            for(size_t k=0;k<childSize();k++)
                child(k)->axpy3(A,X->child(child(k)->_jVec),Y->child(child(k)->_iVec));
        }
        // post-multiply by tensor factor
        else if(_tmp->idx()==iIndex->child(0)){
            for(size_t k=0;k<childSize();k++){
                _tmp->setToZero();
                child(k)->axpy3(A,X->child(child(k)->_jVec),_tmp);
                for(size_t i=0;i<iIndex->childSize();i++){
                    Y->child(i)->axpy((*_tensor)(i,child(k)->_jVec),_tmp);
                }

            }
        }
        // pre-multiply by tensor factor
        else if(_tmp->idx()==jIndex->child(0)){
            for(size_t k=0;k<childSize();k++){
                _tmp->setToZero();
                for(size_t j=0;j<jIndex->childSize();j++){
                    _tmp->axpy((*_tensor)(child(k)->_iVec,j),X->child(j));
                }
                child(k)->axpy3(A,_tmp,Y->child(child(k)->_iVec));
            }
        }
        else DEVABORT(Str("error at")+iIndex->basis()->str()+"<-"+jIndex->basis()->str()
                      +"tmp"+_tmp->idx()->basis());
    }
}

std::string OperatorMap::strNode(int Level) const {
    int digits=Level==Tree_withPtrs?0:Level;
    Str s("","");
    if(iIndex==0)s=s+"F";
    else s=s+_iVec+","+_jVec;
    if(_tensor==0)s=s+" [-]";
    else s=s+" ["+_tensor->rows()+"X"+_tensor->cols()+"] ";
    if(iIndex!=0)s=s+" "+iIndex->axisName()+": "+iIndex->basis()->str()+" <- "+jIndex->basis()->str();
    if(Level==Tree_withPtrs)s+=" ("+tools::str(iIndex)+" <- "+tools::str(jIndex)+")";
    if(_tensor!=0 and Level>0 and Level<15){
        s=s+EigenTools::str(*_tensor,digits);
    }
    return std::move(s);
}
