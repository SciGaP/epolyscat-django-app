// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __COEFFICIENTS__
#define __COEFFICIENTS__

#include "linSpaceVector.h"
#include "tree.h"
#include "labelled.h"
#include <complex>
#include <deque>
#include <set>
#include <unordered_map>
#include "qtEigenDense.h"


//forward declarations
class CoefficientsFloor;
class Discretization;
class Index;
//class Operator;
class Parallel;
/** @defgroup Coefficients
 * \ingroup DiscretizationClasses
 * \brief Tree-structures that hold coefficient data
 *  @{
*/
/// \brief Tree-structered set of expansion coefficients (main class)
class Coefficients : public LinSpaceHilbert<Coefficients>, public Tree<Coefficients>
{
//    friend class Operator;
    friend class Parallel;
    friend class CoefficientsGlobal;
    friend class CoefficientsLocal;
    friend class CoefficientsViewDeep;
    friend class CoefficientsPermute;

    static std::set<const Coefficients*> _viewFloor;
    static std::unordered_map<const Coefficients*,std::vector<std::complex<double> >* > _centralStorage;

    // this could go back into Lablled<...>
    static std::map<const Coefficients*, std::string> _labels;
    static long unsigned int _labelCurrent;

    std::complex<double>* _cData; ///< pointer to data, if !=0 data is Index-ordered, ie. sorted as _cIndex at node
    void setOrderedData(std::complex<double> * CData=0);        ///< set the _cData pointers
    void treeOrderStorage(std::complex<double>*& NextData); ///< rearranged storage sequential in tree (for CoefficientsLocal not Index-ordered)
    std::vector<std::complex<double> >  * nodeStorage() const; ///< (OBSOLESCENT) contiguous storage on present node - created if none, floor points into this

    void setFloorData(std::complex<double> * Data);
    unsigned int storageSize() const;
    void storageAssign(int Size, std::complex<double> Val) const;
    void storageAssign(const std::vector<std::complex<double> > & Stor) const;
    void nodeStorageClear() const;


protected:
    const Index* _cIndex;         ///< contains all structural information for the coefficients
    /// re-arrange storage to tree-contiguous from present node downward
    void clearStorage(); ///< set storage size=0, leave storage pointers (should have been re-assigned before)
    void replaceStorage(std::vector<std::complex<double> > & Stor); ///< replace storage with Stor, leave storage pointers (should have been re-assigned before)
    void averageMargin(Coefficients* Upper, unsigned int Level, double Scale); ///< average margin values on Level
    /// average margin values on Level / optoinally scale values by Scale (needed when dual vectors are expanded)
    void averageFloorMargin(const Index * ILow, std::complex<double>* CLow, const Index * IUpp, std::complex<double>* CUpp, unsigned int Level, const double Scale);
    void makeView() const;// {_viewFloor.insert(this);}
    std::string setToFunctionGrid(std::string Function); ///< evaluate expression on grid valeus for basisIntegrable
    virtual bool notNull() const {return true;}; ///< determines wether true branch
    Coefficients(const Index *I, std::complex<double> Val, std::complex<double> *CDataBegin); ///< initialize to Val
    Coefficients(const Coefficients& Other, std::complex<double> *CDataBegin /** start address of contiguous Coefficient storage*/);
    Coefficients(const Index *I, const Index *IOrigin, unsigned int FloorDepth, std::complex<double>* CDataBegin); ///< continue below the storage to HeightAboveFloor (
public:
    static std::string DEBUGstatus;
    static const int floorPointers=Tree_ptrsOnly-1;

    std::string hash() const;
    bool hasFloorData() const;  ///< non-zero data on floor level
    bool isView() const {return _viewFloor.end()!=_viewFloor.find(this);}
    const Index* idx() const {return _cIndex;} ///< Index tree
    void setIdx(const Index* Idx){_cIndex=Idx;}///< set Index tree for curren Coefficients

    // manipulate storage
    void treeOrderStorage(){std::complex<double>*NextData=0; treeOrderStorage(NextData);}
    // ====== constructors =====================================================
    Coefficients();
    Coefficients(const Coefficients& Other):Coefficients(Other.idx()){*this=Other;}
    Coefficients(const Index *I, std::complex<double> Val=0.); ///< initialize to Val, guarantee contiguous storage
    Coefficients(int FloorDepth, const Index *I,  std::complex<double> Val=0., std::complex<double> *CDataBegin=0); ///< alternate constructor, with explicit floor level
    Coefficients(const Index *I, const std::vector<std::complex<double>> &Vals); ///< initialize to all different vectors
    Coefficients(std::string File, const Index *& NewIndex); ///< read NewIndex and matching Coefficients
    Coefficients(const Index *I, Coefficients * Origin); ///< construct a view on Origin by Index I

    /// view on Origin, place views where SelectView(node) is true
    ///<br>use, e.g., as: Coefficients(Origin, [&](Coefficients* Node){return Node->isLeaf();})
    Coefficients(Coefficients& Origin, std::function<bool(Coefficients*)> SelectView);

    virtual ~Coefficients();
    // ====== functions =====================================================

    // LinSpaceHilbert pure virtual functions
    double norm() const;       ///< norm, NOTE: !!! max(|real|,|imag|)-norm, NOT L2-norm !!!
    double pNorm(int P) const;       ///< p-norm P=0: infty norm, P=1: 1-norm, P=2: [L2-Norm]^2
    virtual long unsigned int size() const;
    virtual Coefficients& axpy (std::complex<double> A, const Coefficients &X, std::complex<double> B);
    virtual Coefficients & operator+=(const Coefficients &rhs);
    virtual Coefficients & operator-=(const Coefficients &rhs);
    virtual std::complex<double> dotProduct(const Coefficients &Other) const {return innerProduct(&Other,false);}
    virtual std::complex<double> scalarProduct(const Coefficients &RightHandVector) const;
    virtual Coefficients & operator*=(std::complex< double > A);

    virtual Coefficients& cwiseMultiply(const Coefficients &B );/// replace this[i] <- this[i]*B[i]
    virtual Coefficients& cwiseDivide(const Coefficients &B);  /// replace this[i] <- this[i]/B[i]
    virtual Coefficients& cwiseRelativeError(const Coefficients &B,double Bound=1.);  /// replace this[i]<-2*(this[i]-B[i])/max(this[i]+B[i],Bound)

    // convenience linear algebra
    virtual Coefficients & operator=(const Coefficients &rhs);
    void axpy(std::complex<double> A, const Coefficients &X){axpy(A,&X);} /// pure axpy
    virtual void axpy(std::complex<double> A, const Coefficients *X); /// pure axpy (with pointer argument)
    virtual void axpy(const Eigen::MatrixXcd & Amat, const Coefficients *X); /// axpy with a matrix (only if sub-vectors are equivalent)

    /// if pseudoScalar==true, do not complex conjugate in inner Product
    virtual std::complex<double> innerProduct(const Coefficients* ket, bool pseudoScalar=false) const;

    std::complex<double> dotProduct(const Coefficients & Ket){return innerProduct(&Ket,false);}

    /// inner product restriced to unscaled region (=0 in scaled region)
    virtual std::complex<double> innerProductUnscaled(const Coefficients* ket) const;
    void scale(const std::complex<double> A);  ///< pure in-place scale

    // vector info
    bool isZero(double eps=0.) const ; ///< true if all |coefficients| <=eps;
    bool isNan() const;
    std::complex<double> cMaxNorm() const; ///< coefficient with maximum norm
    double maxCoeff() const {return std::abs(cMaxNorm());}          ///< L-infty norm: max[i] |xData()[i]|
    std::vector<std::complex<double>> values() const;

    // assign and manipulate coefficients
    void reset(const Index* Idx);    ///< reset Coefficients to Idx structure
    void purgeNearZeros(double Eps); /// set coefficiens<Eps == 0
    void setToZero();                ///< fill coefficient with zeros (optionally set a dummy norm=1
    void setToConstant(std::complex<double> Val); ///< fill coefficient with number
    void setToRandom();              ///< fill coefficient with random numbers (use srand48(int); before!)
    /// fill coefficients with a given function of coordinates (SLOW!)
    ///
    /// Example: c.setToFunction("pow[2](Eta)*exp(-Rn)")
    void setToFunction(std::string Function);
    void conjugate();                ///< conjuagate the coefficient vector

    /// vector of pointers to the actual coefficients; pC.size() = numberCoeffs();
    void pointerToC(std::vector<std::complex<double> *> & pC);

    // continuity conditions, topology
    virtual void makeContinuous(double Scale=1.); ///< impose continuity conditions; optionally Scale margin values (needed when dual vectors are expanded)
    std::complex<double>* data(); ///< storage address (=0 if no storage)
    const std::complex<double>* data() const; ///< storage address (=0 if no storage)

    virtual std::complex<double>* storageData() const;

    // arrange data in address space
    void unsetOrderedData(int Level=INT_MAX); ///< remove orderedData() pointers down to Level
    bool isContiguous() const; ///< true if all data below node is stored contiguously (not necessarily Index-ordered)

    // input/output
    using Tree::str;
    std::string strNode(int Precision=0) const; ///< node data (for Tree)
    void write(std::ofstream & Stream, bool Header, std::string Kind="IndexStructure") const; ///< binary write
    void print(std::ofstream & stream) const;        ///< ascii write
    bool read(std::ifstream & stream, bool header);  ///< binary read
    static void print(const std::vector<Coefficients *> Coeff, std::string Text);

    // Tree virtual functions
    void nodeCopy(const Coefficients* Other, bool View);
    bool nodeEmpty() const;
    bool nodeEquivalent(const Coefficients *Other) const;

    // access to Coefficients data
    std::complex<double> * anyData() const;     ///< pointer to data, not necessarily Index-ordered (=0 if not available)
    std::complex<double> * orderedData() const; ///< pointer to Index-ordered data (=0 if not available)
    std::complex<double> * floorData() const;   ///< pointer to Index-ordered data in floor (=0 if not floor)
    std::complex< double > floorInnerProduct(const Coefficients *ket, bool pseudoScalar) const;

    // tree-related functions
    Coefficients *firstFloor() const ; ///< first leaf below present Coefficient;
    Coefficients* retrieve(const Index* Idx);///< coefficient at index

    void examplePermute(); ///< demonstrates how to permute
    static void cleanUp(); ///< for final cleanup, shuts up valgrind about _centralStorage data

    void plot(std::string PltFile="");
};

// left hand scalar functions
inline Coefficients operator*(std::complex<double> A, const Coefficients & C){Coefficients Ca(C); Ca*=A; return Ca;}

/** @} */ // end group Coefficients
#endif
