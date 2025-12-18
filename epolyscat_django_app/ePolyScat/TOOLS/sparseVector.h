// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef SPARSEVECTOR_H
#define SPARSEVECTOR_H

#include <vector>
#include "tree.h"

template<class T>
class Scalar{
public:
    typedef T scalar;    ///< T type
    T data; ///< scalar data
public:
    Scalar():data(0){}
    Scalar(const T &Data):data(Data){}
    std::string str() const {return tools::str(data);}
    T & val(unsigned int I){ return data;}
    void axpy(T A, T & X){data+=A*X.data;}
    void scal(T A){data*=A;}

    unsigned int size() const {return 1;}
    double normSquare() const {return std::abs(data*data);}
};

template<class T>
class Vector{
public:
    typedef T scalar;    ///< T type
    std::vector<T> data; ///< vector data
public:

    Vector(){}
    Vector(const std::vector<T>&Data):data(Data){}

    std::string str() const {return tools::str(data,",",4);}
    T & val(unsigned int I){ return data[I];}

    void scal(T A){
        if(A==0 or data.size()==0)return;
        for(typename std::vector<T>::iterator y=data.begin();y!=data.end();y++)*y*=A;
    }

    void axpy(T A, const Vector<T> & X){
        if(A==0. or X.data.size()==0)return;
        typename std::vector<T>::const_iterator x=X.data.begin();
        if(data.size()==0){
            data=X.data;
            scal(A);
        } else {
            for(typename std::vector<T>::iterator y=data.begin();y!=data.end();y++,x++)*y+=*x*A;
        }
    }
    unsigned int size() const {return data.size();}
    double normSquare() const {
        double normSq=0.;
        for(typename std::vector<T>::const_iterator y=data.begin();y!=data.end();y++)normSq+=std::abs((*y)*(*y));
        return normSq;
    }
};

/// \ingroup Linalg
///@brief sparse vector is a tree of indices and data
///
/// prefer use where tree overhead is negligible
///
/// the index is unique for all children of a given parent node
/// data must allow vector space operations and provide an L2-norm
/// data item for SparseVector: scalar entry
template <class T>
class SparseVector:public IndexTree<SparseVector<T> >
{
    // let any try and all fellow SparseVector access the data
    template <class U> friend class Tree;
    template <class V> friend class SparseVector;

    //    int idx; ///< index of the node, unique among children of the same parent
    bool recalc;         ///< data changed, need to recalculate auxiliary data

    double normSq;       ///< L2-norm of tensor
    unsigned int treeSz; ///< number of non-zero elements in tree

    T data;
protected:
    static SparseVector<T> pathIndex(const std::vector<int> & IndexPath); ///< single path with indices as in IndexTree (for insert, erase)
    std::string strNode() const;

public:
    SparseVector<T>(int Index=0):IndexTree<SparseVector<T> >(Index),recalc(true){}
    SparseVector<T>(int Index, const T & Data):IndexTree<SparseVector<T> >(Index),recalc(true),data(Data){}

    typename T::scalar val(unsigned int I){return data.val(I);} ///< return I'th scalar data from present node

    void matrix(int I, int J, const typename T::scalar & Data);///< insert matrix element

    SparseVector<T> &scal(typename T::scalar A); ///< scale by A

    /// vector axpy operation: this+=A*X
    SparseVector<T> &axpy(typename T::scalar A, const SparseVector<T>& X);

    /// matrix multiplication onto 1st index; skip contributions with L2-norm < Eps
    SparseVector<T> &addMatrixVector(const SparseVector<Scalar<typename T::scalar> > &Mat, const SparseVector<T> &Vec, double Eps=0.);

    void cwiseProduct(const SparseVector<T> & X); ///< component-wise product of tensors

    double normSquare() const     {if(not recalc)return normSq; ABORT("tensor data outdated, run properties()");} ///< squared L2-norm
    unsigned int dataSize() const {return data.size();} ///< data entries on node
    unsigned int treeSize() const {if(not recalc)return treeSz; ABORT("tensor data outdated, run properties()");} ///< data in and below node
    SparseVector<T> & properties(); ///< recalculate properties stored with tensor

    SparseVector<T> & purge(double Eps); ///< remove near-zeros such that total normSq changes by at most Eps
    SparseVector<T> & zero(const std::vector<int> & IndexPath=std::vector<int>()); ///< set subtree at index position =0 (=remove from tree), return false if not in tree

    static void Test();
};

template <class T>
SparseVector<T>& SparseVector<T>::scal(typename T::scalar A){
    if(A==1.)return *this;
    recalc=true;
    data.scal(A);
    // for reasons that I do not care to understand (rather I hope for a better language than C++)
    // here, we cannot bluntly refer to menbers by their name, but need the "this"
    for(unsigned int k=0;k<this->childSize();k++)this->child(k)->scal(A);
    return *this;
}

template <class T>
SparseVector<T>& SparseVector<T>::axpy(typename T::scalar A, const SparseVector<T> & X){
    if(A==0.)return *this;
    recalc=true;
    data.axpy(A,X.data);
    typename std::vector<SparseVector<T>*>::iterator p=this->child.begin();
    unsigned int ip=0;
    for(unsigned int k=0;k<X.childSize();k++){
        if(not this->indexPos(X.child(k)->index(),ip)){
            this->attach(new SparseVector(*X.child(k)));
            this->child(ip)->scal(A);
        } else {
            this->child(k)->axpy(A,*X.child(k));
        }
    }
    return *this;
}

template <class T>
SparseVector<T> & SparseVector<T>::addMatrixVector
(const SparseVector<Scalar<typename T::scalar> > & Mat, const SparseVector<T> &X, double Eps){

    // CAUTION: this is not using the matrix norm!
    if(Mat.normSquare()*X.normSquare()<=Eps)return *this; // skip near-zero contributions
    recalc=true;

    if(this==&X)ABORT("no aliasing allowed here");

    typename std::vector<SparseVector<Scalar<typename T::scalar> >*>::const_iterator iMat,ijMat;
    //    typename std::vector<SparseVector<T>*>::iterator y=this->child.begin();
    unsigned int iy=0;
    iMat=Mat.child.begin();
    for(iMat=Mat.child.begin();iMat!=Mat.child.end();iMat++){
        if((**iMat).normSquare()*X.normSquare()>Eps){

            if(not this->indexPos((**iMat).index(),iy))
                this->attach(new SparseVector<T>((**iMat).index()));

            unsigned int ix=0;
            for(ijMat=(**iMat).child.begin();ijMat!=(**iMat).child.end();ijMat++)
                if(X.indexPos((**ijMat).index(),ix))
                    this->child(iy)->axpy((**ijMat).val(0),*X.child(ix)); // axpy of children
        }
    }

}

template <class T>
void SparseVector<T>::matrix(int I, int J, const typename T::scalar & Data){
    this->attach(new SparseVector<T>(J, Scalar<typename T::scalar>(Data)),
                 SparseVector<T>::pathIndex(std::vector<int>(1,I)));}

template <class T>
std::string SparseVector<T>::strNode() const {
    std::string s=" "+tools::str(this->index())+"(";

    if(recalc)s+="*|*";
    else s+=tools::str(this->treeSize())+"|"+tools::str(this->normSquare(),3);
    s+=")";
    if(this->data.size()>0)s+=": "+data.str();
    return s;
}

template <class T>
SparseVector<T> & SparseVector<T>::properties(){
    if(not recalc)return *this;
    treeSz=this->data.size();
    normSq=this->data.normSquare();
    for(typename std::vector<SparseVector<T>*>::const_iterator x=this->child.begin();x!=this->child.end();x++){
        (**x).properties();
        treeSz+=(**x).treeSize();
        normSq+=(**x).normSquare();
    }
    recalc=false;
    return *this;
}

template <class T>
SparseVector<T> SparseVector<T>::pathIndex(const std::vector<int> &IndexPath){
    SparseVector<T> p,*pp;
    pp=&p;
    for(unsigned int k=0;k<IndexPath.size();k++){
        pp->child.push_back(new SparseVector<T>(IndexPath[k]));
        pp=pp->child(0);
    }
    std::cout<<"path "<<p.child(0)->childSize()<<"\n"<<p.str()<<std::endl;
    return p;
}

#endif // SPARSEVECTOR_H
