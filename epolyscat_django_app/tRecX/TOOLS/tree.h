// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef TREE_H
#define TREE_H

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "tools.h"
#include "str.h"

// to be done:
// add lots of checks for in debug mode

#include "log.h"

const int Tree_noNodeInfo=-INT_MAX;
const int Tree_defaultKind=-1;
const int Tree_ptrsOnly=-INT_MAX+1;
const int Tree_ptrsSizes=-INT_MAX+200;
const int Tree_withPtrs=Tree_defaultKind+100;

///@brief a lean Tree class (thanks to Andreas Swoboda for the "Curiously Recurring Template Pattern" design)
///
/// Usage example in class SparseVector
template <typename T>
class Tree
{
    template <class U> friend class Tree;
    template <class U> friend class TreeDerived;


    static std::string indent;

    /// increment Pos in child while it is below Node, true if Pos at same position as Node
    bool childPos(const T* Node, typename std::vector<T*>::iterator &Pos) const;

    // auxiliary routine for inserting or appending an element to child
    inline typename std::vector<T*>::iterator insertChild(typename std::vector<T*>::iterator Pos, T* NewChild); // insert at Pos among child


    const T* _parent;
    std::vector<T*> * _child;

    static std::set<const T*> _keepChild;
    bool keepChild() const {return _keepChild.end()!=_keepChild.find(const_cast<T*>(dynamic_cast<const T*>(this)));}

    typename std::vector<T*>::iterator childEnd()  const {return _child->end();}
    typename std::vector<T*>::iterator childBegin()const {return _child->begin();}

    void childPushBack(T* Child){if(_child==0)_child=new std::vector<T*>; _child->push_back(Child); } ///< append sub-node


protected:


    // use of virtual functions:
    // in the template, call of virtualFunction(...) or this->virtualFuction(...) will call the base version,
    // WARNING: advice against such use, because it rarely really is the intended effect
    virtual std::string strNode(int Level=Tree_ptrsOnly) const; ///< string showing data of node

public:

    long diagnoseNodeCount() const; ///< total number of nodes in Tree
    size_t diagnoseSizeOf() const;
    virtual size_t diagnoseSizeOfNode() const;

    Tree(const T* Parent=0):_parent(Parent),_child(0){} ///< elementary constructor
    Tree(const Tree<T> &Other, bool View=false); ///< (deep) copy constructor
    Tree<T> & operator=(const Tree<T> & Other); ///< (deep) assignement operator

    T* deepCopy(size_t PruneDepth=INT_MAX) const; ///< return a deep copy down to depth()<=PruneDepth
    T* view() const; ///< return new tree with nodeCopy(this,view) and view on children

    Tree<T>(const Tree<T> * InTree, std::vector<unsigned int> Permute); ///< return a permuted view of a tree

    /// return a tree with levels permuted
    T & permute(std::vector<unsigned int> Perm /** permutation of levels */, T & Out, bool View=false /** =true: copy only permuted nodes */ ) const;
    /// permute levels in tree
    void permuteInPlace(std::vector<unsigned int> Perm /** permutation of levels */);
    T* permute(std::vector<unsigned int> Perm /** permutation of levels */,bool View=false) const;
    virtual void purge(unsigned int Height); ///< remove branches that do not have at least one node on Level below present (virtual as some trees need special care)

    T * subTree(bool (*Select)(const T*), T *Out, bool View=false) const; ///< create a copy (or view: View=true) of a tree by criterion Select

    /// new Tree from Levels of this
    ///
    /// aborts if levels to be omitted do not consist of equivalent subtrees
    T * factor(const std::vector<int> Levels /** include Levels */,
               bool Exclude=false /** revert meaning: exclude Levels */) const;

    /// true if structure matches and all nodes are nodeEquivalent
    bool treeEquivalent(const T * Other, const std::vector<int> OnlyLevels) const;

    /// true if all branches are treeEquivalent
    bool branchesEquivalent() const;


    // various relatives
    T * root() const {if(parent()==0)return const_cast<T*>(dynamic_cast<const T*>(this));return parent()->root();} ///< first node in tree

    virtual const T * parent() const {return _parent;}          ///< parent node
    virtual const T *&parentRef() {return _parent;}          ///< reference to parent

    T * child(unsigned int N) const {return _child->data()[N];} ///< return N'th child
    T*& child(unsigned int N) {return _child->data()[N];} ///< return N'th child
    T*& childRef(unsigned int N) {return _child->data()[N];} ///< return reference to N'th child
//    std::vector<T*> const children() {if(not _child)return {}; return *_child;}
    std::vector<T*> const children() const {if(not _child)return {}; return *_child;}
    unsigned int childSize() const {if(_child==0) return 0; return _child->size();} ///< number of sub-nodes
    T * childBack() const {return _child->back();} ///< return last child
    void childPop() {delete _child->back(),_child->pop_back();}
    void childAdd(T* Child){childPushBack(Child); _child->back()->parentRef()=dynamic_cast<const T*>(this);} ///< append sub-node
    void childInsert(size_t Pos, T* Child){if(not _child)_child=new std::vector<T*>();insertChild(childBegin()+Pos,Child);} ///< append sub-node
    void childView(T* Child){childPushBack(Child);_keepChild.insert(_child->back());} ///< append sub-node
    void childReplace(unsigned int N, T* Child){delete child(N);child(N)=Child;child(N)->parentRef()=dynamic_cast<const T*>(this);} ///< replace child(N) with Child, data of previous child(N) is deleted
    void childErase(unsigned int N){delete child(N);childEraseNode(N);} ///< delete child, removed child entry
    void childEraseNode(unsigned int N){_child->erase(childBegin()+N);if(childSize()==0){delete _child; _child=0;}} ///< remove node, do not delete subtrees
    void childDetach(unsigned int N, T* &Detached){Detached=child(N);childEraseNode(N);} ///< return reference to child, remove reference from parent tree
    void childrenMove(T &NewParent); ///< move all children to new parent

    T * rightSibling() const; ///< return right sibling, =0 if none
    T * nodeRight(const T *Root=0) const;    ///< node to the right of this in subtree starting from Root (Root==0: complete tree); return 0 if last node on level; error if this not desendent of Root
    T * nodeNext(const T* Root=0) const;    ///< next right, or else, first non next lower level relative to Branch (Root==0: complete Tree); return 0 if last node on level; abort with error if this is not desendent of Branch
    /// leftmost node in Branch, Decscend levels down from present, Branch==0: only in present subtree
    T * descend(unsigned int Descend=1, const T *Branch=0) const;
    /// return ancesctor up by Ascend levels
    T * ascend(size_t Ascend=1) const {const T* a=dynamic_cast<const T*>(this);for(size_t k=0;k<Ascend;k++)a=a->parent();return const_cast<T*>(a);}; /// leftmost node in Branch, Decscend levels down from present, Branch==0: only in present subtree
    /// remove node from parent, return with parent set to 0,
    T * detach() const {T* res=const_cast<T*>(dynamic_cast<const T*>(this));
                        if(parent())const_cast<T*>(parent())->childEraseNode(nSibling());
                        return res;}

    T* nodeAt(const std::vector<unsigned int> & Idx) const {
        if(Idx.size()==0)return const_cast<T*>(dynamic_cast<const T*>(this));
        if(Idx[0]>=childSize())return 0;
        if(Idx.size()==1)return child(Idx[0]);
        return child(Idx[0])->nodeAt(std::vector<unsigned int>(Idx.begin()+1,Idx.end()));
    }

    int nSibling() const; ///< position of node among children of its parent()
    std::vector<unsigned int> index() const; ///< sibling numbers from top to including present

    unsigned int size(unsigned int Depth) const; ///< number of nodes at Depth
    unsigned int levelSize() const; ///< number of nodes on level
    unsigned int leafSize() const;  ///< number of leafs in tree
    unsigned int depth() const;     ///< number of links to root
    unsigned int height() const;    ///< number of links firstLeaf()
    unsigned int levelRank(const T* Parent=0) const; ///< number of nodes to the left, on present level

    bool isLeaf() const; ///< a Leaf is a tree without children
    bool isSubtree(const T* Root) const; ///< true if this is subTree of Root
    T *firstLeaf() const; ///< leaf at the left edge of tree
    T *nextLeaf(const T * Subtree=0) const; ///< next leaf (to the right) at lower edge of the tree starting from Subtree

    // for rearranging trees, class T must tell which data it considers relevant and when it is to be considered empty
    // functions below all refer to the given node, not to its relations in the tree
    virtual bool nodeEmpty() const {DEVABORT("need to define nodeEmpty");};
    virtual bool nodeEquivalent(const T* Other) const {DEVABORT("need to define nodeEquivalent");return Other==0;}; // shut up unused variable warnings
    virtual void nodeCopy(const T* Node, bool View) {DEVABORT("need to define nodeCopy");if(View)_child=Node->_child;}; // shut up unused variable warnings

    typename std::vector<T*>::iterator attach(T *Pointer,  const T & insertPath=T()); ///< attach Tree POINTER after the end of insertPath

    bool isView() const {return keepChild();}

    virtual bool strDataShow(){return true;}
    virtual std::string str(int Kind=Tree_defaultKind /** see strNode(...) re-implementation for meaning (default: node pointer) */,
                            int Depth=INT_MAX /** print for depth()<=Depth */ ) const;  ///< default printing

    /// remove all sub-trees
    void clear();
    
    /// All children at a given depth below current node ordered, be sure to give an empty vector!
    void childrenAtDepth(int depth, std::vector<const T*>& children) const;

    std::vector<int> sizeSummary() const;

public:
    // virtual destructor makes sure that before this destructor is called, the destructor of the derived class is called
    // virtual is inheritable, i.e. destruction starts from the last derived class and works its way back
    virtual ~Tree();
};

///@brief a Tree where each node carries an integer index
template<class T>
class IndexTree:public Tree<T>{
    template<class U> friend class Tree;

protected:
    // these should be used by Tree only:
    inline bool matchChild(typename std::vector<T*>::iterator &Pos) const {return idx==(**Pos).idx;}
    inline bool aboveChild(typename std::vector<T*>::iterator &Pos) const {return idx >(**Pos).idx;}

    const int idx; // do not manipulate index after construction
    inline int index() const {return idx;}

    /// advance directly for index matching (no detour through matchChild/aboveChild
    bool indexPos(int Index, unsigned int & Pos) const  {
        while( Pos < this->childSize() and Index >this->child(Pos)->idx)Pos++;
        return Pos < this->childSize() and Index==this->child(Pos)->idx;
    }

public:
    IndexTree(int Index=0):idx(Index){}
};

//=========================================================================================
//--- codes -------------------------------------------------------------------------------
//=========================================================================================

template <class T>
std::set<const T*> Tree<T>::_keepChild;

template <class T>
Tree<T>::~Tree(){
    bool del=not keepChild();
    for(unsigned int k=0;del and k<childSize();k++)delete child(k);
    delete _child;
    _child=0;
    if(del)_keepChild.erase(dynamic_cast<T*>(this)); // remove from set to be kept
}

template <class T>
void Tree<T>::clear(){
    bool del=not keepChild();
    for(unsigned int k=0;k<childSize();k++)if(del)delete child(k);
    if(_child!=0){
        _child->clear();
        delete _child;
        _child=0;
    }
}

template <class T>
std::string Tree<T>::indent="  ";

template <class T>
Tree<T>::Tree(const Tree<T> &Other, bool View):_child(0){ ///< (deep or view) copy constructor
    if(View)_keepChild.insert(this);
    nodeCopy(dynamic_cast<T*>(this),View);
    for (unsigned int i=0; i!=Other.childSize(); ++i){
        if(View) childAdd(Other.child(i));
        else childAdd(Other.child(i)->deepCopy());
    }
}

template <class T>
void Tree<T>::purge(unsigned int Height){
    if(Height==1)return;
    for(unsigned int k=childSize();k>0;k--){
        //        if(child(k-1)->descend(Height-1)==0)childErase(k-1);
        //        else {
        child(k-1)->purge(Height-1);
        //        }
    }
    for(unsigned int k=childSize();k>0;k--)
        if(child(k-1)->isLeaf())childErase(k-1);
}


// this should be replaced with deepAssign(..,INT_MAX)
template <class T>
Tree<T> & Tree<T>::operator=(const Tree<T> &Other){ ///< (deep) assignment
    if(this==&Other)return *this;
    _parent=0;
    for (unsigned int i=0; i!=Other.childSize(); ++i) {
        typename std::vector<T*>::iterator pos=childBegin();
        if(childPos(Other.child(i),pos)){
            // overwrite
            delete *pos;
            *pos=new T(*Other.child(i));
        } else {
            // insert/append new
            childPushBack(new T(*Other.child(i)));
        }
    }
}

template <class T>
bool Tree<T>::childPos(const T* Node,typename std::vector<T*>::iterator & Pos) const{
    while( Pos < childEnd() and Node->aboveChild(Pos)) Pos++;
    return Pos < childEnd() and Node->matchChild(Pos);
}
template <class T>
std::string Tree<T>::strNode(int Level) const {
    std::string s;
    if(parent()==0)s="R";
    if(childSize()==0)s="l";
    else s+="+"+tools::str(childSize()) ;
    if(Level==Tree_ptrsOnly)s+=" "+tools::str(this);
    return s;
}

template <class T>
bool Tree<T>::isLeaf() const {return childSize()==0;}

template<class T>
unsigned int Tree<T>::depth() const{
    if(parent()==0)return 0;
    return parent()->depth()+1;
}

template<class T>
unsigned int Tree<T>::height() const{
    if(isLeaf())return 0;
    return child(0)->height()+1;
}


template<class T>
unsigned int Tree<T>::levelRank(const T* Parent) const {
    if(parent()==0 or this==Parent)return 0;

    unsigned int s=nSibling();
    unsigned int d=1;
    if(parent() == Parent) return s;
    const T* p=parent();
    for(; p!=0 and p->parent()!=Parent;p=p->parent()){
        for(int k=0;k<p->nSibling();k++){
            s+=p->parent()->child(k)->Tree::size(d);
        }
        d++;
    }
    if(p==0)ABORT("argument is not parent of node");
    return s;
}

template<class T>
T * Tree<T>::rightSibling() const {
    if(parent()==0)return 0;
    unsigned int n=nSibling();
    if(n==parent()->childSize()-1)return 0;
    return parent()->child(n+1);
}

template<class T>
T* Tree<T>::firstLeaf() const {
    if(isLeaf())return dynamic_cast<T*>(const_cast<Tree<T>*>(this));
    else return child(0)->firstLeaf();
}
template<class T>
T * Tree<T>::nextLeaf(const T* Subtree) const {
    if(this==Subtree)return 0;
    if(not isLeaf())ABORT("nextLeaf can only be called by leaf");

    const T * next=rightSibling();
    if(next!=0){
        next=next->firstLeaf();
    }
    else {
        // ascend until new subtree below Root
        next=parent();
        while(next!=Subtree and next->rightSibling()==0){
            next=next->parent();
        }
        if(next==Subtree)next=0;
        if(next!=0)next=next->rightSibling()->firstLeaf();
    }
    return const_cast<T*>(next);
}

template<class T>
unsigned int Tree<T>::levelSize() const {
    return root()->Tree::size(depth());
}

template<class T>
unsigned int Tree<T>::size(unsigned int Depth) const {
    if(Depth<2){
        if(Depth==0)return 1;
        return childSize();
    }
    unsigned int s=0;
    for(int k=0;k<childSize();k++)s+=child(k)->Tree::size(Depth-1);
    return s;
}

template<class T>
T * Tree<T>::nodeNext(const T* Branch) const {
    const T* branch=Branch;
    if(branch==0)branch=root();
    T* next=nodeRight(branch);
    if(next==0)
        next=branch->descend(depth()-branch->depth()+1,branch);
    return next;
}

template<class T>
T * Tree<T>::descend(unsigned int Descend, const T* Branch) const {

    if(Branch==0)Branch=dynamic_cast<const T*>(this);

    if(Descend==0)return dynamic_cast<T*>(const_cast<Tree<T>*>(this));
    if(childSize()>0)return child(0)->descend(Descend-1,Branch);
    if(Branch==0)return 0; // no branch specified - left edge only

    T* right=nodeRight(Branch);
    if(right!=0 and not right->isLeaf())return right->descend(Descend,Branch);
    return 0;
}

template<class T>
bool Tree<T>::isSubtree(const T* Root) const {
    if(Root==this)return true;
    const T * s=dynamic_cast<const T*>(this);
    while(s!=0 and s!=Root){
        s=s->parent();
    }
    return s!=0;
}

template<class T>
T * Tree<T>::nodeRight(const T* Root) const {
    if(Root==0)Root=root();
    if(parent()==0)return 0;
    if(this==Root) return 0; // at root of subtree

    T* right=rightSibling();
    if(right!=0)return right;

    unsigned int lev=0;
    const Tree<T>* up=parent();
    while(up->parent()!=0 and up->parent()->_child->back()==up){
        if(up==Root)return 0; // cannot move further up
        up=up->parent();
        if(up==0)ABORT("index "+tools::str(this)+" not in subtree: \n"+Root->Tree::str());
        lev++;
    }
    if(up==Root or up->rightSibling()==0)return 0; //
    return up->rightSibling()->descend(lev+1,Root);
}


template <class T>
int Tree<T>::nSibling() const {
    int nsib=0;
    if(parent())nsib=std::find(parent()->childBegin(),parent()->childEnd(),this)-parent()->childBegin();
    return nsib;
}

template <class T>
std::vector<unsigned int> Tree<T>::index() const {
    if(parent()==0)return std::vector<unsigned int>(0);
    std::vector<unsigned int> i(parent()->index());
    i.push_back(nSibling());
    return i;
}

// auxiliary routine for inserting or appending an element to child
template <class T>
typename std::vector<T*>::iterator Tree<T>::insertChild(typename std::vector<T*>::iterator Pos, T* NewChild){
    typename std::vector<T*>::iterator pos;
    if(Pos==childEnd()){
        childPushBack(NewChild);
        pos=childEnd();
        pos--;
    } else {
        if(_child==0)_child=new std::vector<T*>();
        pos=_child->insert(Pos,NewChild);
    }
    (**pos)._parent=dynamic_cast<const T*>(this);
    return pos;
}

template <class T>
typename std::vector<T*>::iterator Tree<T>::attach(T *Pointer, const T &insertPath)
{
    typename std::vector<T*>::iterator pos=childBegin();
    if(insertPath.isLeaf()){
        if(childPos(Pointer,pos))ABORT("cannot insert, position occupied:"+Pointer->str()+"\nin tree \n"+str());
        pos=insertChild(pos,Pointer);
    } else {
        // create new node and descend
        if(not childPos(insertPath.child(0),pos))pos=insertChild(pos,new T(*insertPath.child(0)));
        (**pos).attach(Pointer,*insertPath.child(0));
    }
    return pos;
}

template <class T>
std::string Tree<T>::str(int Precision, int Depth) const
{
    std::string s=parent()?tools::str(index(),3,"  ")+" |":"|";
    if(isView())s+="V";
    s+="....";

    if(Precision==Tree_ptrsOnly)s+=tools::str(this);
    if(Precision==Tree_withPtrs or Precision==Tree_ptrsSizes)s+=tools::str(this)+": "+strNode(Precision);
    else                        s+=strNode(Precision);

    for(unsigned int k=0;k<childSize() and Depth>0;k++)s+="\n"+child(k)->Tree<T>::str(Precision,Depth-1);
    return s;
}

/// build a copy (or View) of the tree where Select(node) is true
template <class T>
T* Tree<T>::subTree(bool (*Select)(const T *), T *Out, bool View) const {
    if(View)_keepChild.insert(Out);
    Out->nodeCopy(dynamic_cast<const T*>(this),View);
    for(unsigned int k=0;k<childSize();k++){
        if(Select(child(k))){
            Out->childAdd(new T());
            child(k)->subTree(Select,Out->childBack(),View);
        }
    }
    return Out;
}

template <class T>
T* Tree<T>::factor(const std::vector<int> Levels, bool Exclude) const{

    // skip node: (not in list XOR Exclude) and not end of Tree
    if(not isLeaf() and ((std::find(Levels.begin(),Levels.end(),depth())==Levels.end())!=Exclude)){
        if(not branchesEquivalent())std::cout<<"WARNING: omitting node with non-equivalent branches: "+strNode(0)<<std::endl;
        return child(0)->factor(Levels,Exclude);
    }

    // copy node and continue factoring branches
    T* idx=new T();
    dynamic_cast<Tree<T>*>(idx)->nodeCopy(dynamic_cast<const T*>(this),false);
    for(int k=0;k<childSize();k++)idx->childAdd(child(k)->factor(Levels,Exclude));
    return idx;

}

template <class T>
bool Tree<T>::treeEquivalent(const T * Other, const std::vector<int> OnlyLevels) const{

    //    fail=Str("indexSize")+childSize()+Other->childSize();
    if(childSize()!=Other->childSize())return false;

    //    fail="not node equivalent";
    if((OnlyLevels.size()==0 or std::find(OnlyLevels.begin(),OnlyLevels.end(),depth())!=OnlyLevels.end())
            and not nodeEquivalent(Other))return false;

    // equivalence of sub-indices
    //    fail="subIdx";
    for(unsigned int k=0;k<childSize();k++)
        if(not child(k)->treeEquivalent(Other->child(k),OnlyLevels))return false;
    return true;
}

template <class T>
bool Tree<T>::branchesEquivalent() const {
    for (unsigned int k=1;k<childSize();k++)
        if(not child(0)->treeEquivalent(child(k))){
            return false;
        }
    return true;
}


template <class T>
void Tree<T>::permuteInPlace(std::vector<unsigned int> Permute) {
    T *t=this->permute(Permute,false);
    nodeCopy(t,false);
    while(childSize())childPop();
    Tree<T> *tt=dynamic_cast<Tree<T>*>(t);
    tt->childrenMove(*dynamic_cast<T*>(this));
    delete tt;
}

template <class T>
T* Tree<T>::permute(std::vector<unsigned int> Permute,bool View) const {
    T* t=new T();
    this->permute(Permute,*t,View);
    return t;
}


template <class T>
T& Tree<T>::permute(std::vector<unsigned int> Permute, T& Out, bool View) const{

    // check whether Permute is actually a permutation
    for(int k=0;k<Permute.size();k++)
        if(std::count(Permute.begin(),Permute.end(),k)!=1)
            DEVABORT(Str("not a permutation of 0...","")+(Permute.size()-1)+": "+Permute);

    // remove trivial permutations
    while(Permute.size()>0 and Permute.back()==Permute.size()-1)Permute.pop_back();

    T * l=descend(Permute.size());
    if(l==0)ABORT("tree height smaller than permutation size");

    // descend and create path to childView
    T* res=&Out;
    while (l!=0){
        T* kLevel=res;
        std::vector<unsigned int> idx(l->index());

        for(unsigned int k=0;k<Permute.size();k++){
            for(unsigned int m=kLevel->childSize();m<idx[Permute[k]]+1;m++)kLevel->childAdd(new T());
            // get node from the parent tree level corresponding to new
            T* rep=nodeAt(std::vector<unsigned int>(idx.begin(),idx.begin()+Permute[k]));

            // if node is empty, copy data from the corresponding parent node
            if(kLevel->nodeEmpty())kLevel->nodeCopy(rep,View);

            // else check data equivalence
            else if(not rep->nodeEquivalent(kLevel))
                ABORT(rep->strNode(0)+" != "+kLevel->strNode(0)+"\nCannot rearrange level - node data not equivalent: ");

            kLevel=kLevel->child(idx[Permute[k]]);
        }
        if(View)_keepChild.insert(kLevel);
        if(kLevel->childSize()!=0)ABORT("tree position occupied - cannot place new node");
        kLevel->nodeCopy(l,View);
        for(unsigned int k=0;k<l->childSize();k++){
            if(View)kLevel->childView(l->child(k));
            else    kLevel->childAdd(l->child(k)->deepCopy());
        }
        l=l->nodeRight();
    }
    if(Permute.size()>0)
        res->purge(Permute.size()-1); //< last Level (=Permute.size()-1) must be non-empty
    return Out;
}

template <class T>
void Tree<T>::childrenMove(T &NewParent) {
    if(not NewParent._child)NewParent._child=new std::vector<T*>();
    NewParent._child->insert(NewParent._child->begin(),_child->begin(),_child->end());
    for(T* c: *_child)c->_parent=&NewParent;
    delete _child;
    _child=0;
}

template <class T>
T* Tree<T>::deepCopy(size_t PruneDepth) const {
    T* t=new T();
    t->nodeCopy(dynamic_cast<const T*>(this),false);
    if(not PruneDepth)return t;

    for(unsigned int k=0;k<childSize();k++)
        t->childAdd((child(k)->deepCopy(PruneDepth-1)));
    return t;
}

template <class T>
T* Tree<T>::view() const{
    T* t=new T();
    
    t->nodeCopy(dynamic_cast<T*>(this),true);
    _keepChild.insert(this);
    for(unsigned int k=0;k<childSize();k++)t->childAdd(child(k));
    return t;
}

template<class T>
void Tree<T>::childrenAtDepth(int depth, std::vector<const T*>& children) const{
    if(depth==0){
        children.push_back(dynamic_cast<const T*>(this));
    }else{
        for(unsigned int i=0; i<childSize();i++){
            child(i)->childrenAtDepth(depth-1, children);
        }
    }
}

template<class T>
unsigned int Tree<T>::leafSize() const{
    if(childSize()==0) return 1;
    int res=0;
    for(int i=0; i<childSize(); i++){
        res +=child(i)->leafSize();
    }
    return res;
}

template<class T>
std::vector<int> Tree<T>::sizeSummary() const{
    std::vector<int> result(1, 1);
    for(unsigned int i=0; i<childSize(); i++){
        std::vector<int> tmp = child(i)->sizeSummary();
        for(unsigned int j=0; j<tmp.size(); j++){
            if(result.size()<=j+1) result.push_back(0);
            result[j+1] += tmp[j];
        }
    }
    return result;
}

template<class T>
size_t Tree<T>::diagnoseSizeOf() const{
    size_t siz=diagnoseSizeOfNode();
    for(int k=0;k<childSize();k++)siz+=child(k)->diagnoseSizeOf();
    return siz;
}
template<class T>
size_t Tree<T>::diagnoseSizeOfNode() const{ return sizeof *this;}

template<class T>
long Tree<T>::diagnoseNodeCount() const{
    long cnt=1;
    if (_child!=0)
        for (auto c: *_child)cnt+=c->diagnoseNodeCount();
    return cnt;
}

#endif

