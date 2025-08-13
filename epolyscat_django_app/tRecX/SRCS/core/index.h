// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __INDEX__
#define __INDEX__

#include <exception>

#include "tree.h"
#include "labelled.h"
#include "discretizationDerived.h"
#include "operatorTree.h"
#include "indexOverlap.h"
#include <cstdint>

//forward declaration
class Discretization;
class IndexFloor;
class OperatorSingle;
class UseMatrix;
class Inverse;
class BasisAbstract;
class BasisNdim;
class BasisGrid;
class BasisIntegrable;
class OperatorOverlap;
class OperatorTree;
class IndexConstraint;

//class BasisSet;

/// \defgroup Index
/// \ingroup DiscretizationClasses
/// \brief Recursive indices and derived forms

/// \ingroup Index
/// \brief Index tree (central class of the code)
///
/// Structure:
/// - every Index is the node of a tree
/// - every node has a basis() and an axisName()
/// - usually childSize()==basis()->size()
/// - hasFloor(): from here onward coefficients are guaranteed to be stored contiguously in index order
/// - isBottom(): lowest non-trivial node, childSize()==0 and basis()->size()>0;  (legacy: has childSize()==basis()->size() dummy children)
/// - bottomExpand(): expand bottom to have basis() children (if needed), i.e. bring to legacy form
/// - bottomUnexpand(): remove children from bottom, if they were previously added by bottomExpand()
/// - any constructor should execute sizeCompute()
/// - after sizedCompute(), size() indicates Coefficients::size() starting from a given node

class Index: public Tree<Index>, public Labelled<Index> {
    friend class Discretization;
    friend class Coefficients;
    friend class Gauge;
    friend class IndexOverlap;

    static std::vector<const BasisAbstract*> _allBases;
    static std::vector<std::string> _allAxisNames;
    static std::map<std::string,std::string> _axisSubset;

    Index(std::istream &Stream, std::map<uint8_t,uint8_t> AxisNumber, std::map<uint16_t,uint16_t> BasisNumber,bool Rewind);
    const BasisAbstract * basisFromGrid(int ContractFactor, std::vector<const Index*> Path) const;

    Index* toIndexBasis(std::vector<int> ContractFactor, std::vector<const Index*> Path) const;

protected:
    // here we economize on storage as much as possible
    mutable uint32_t _size;
    uint16_t _indexBas;
    uint8_t _indexAx;
    char _indexKind;

    void contractedNumbering(std::vector<unsigned int> & Numb, std::vector<unsigned int> & Mult,  const Index *From, int &CurMax, int &CurPos) const;

public:
    static bool noDum;
    static std::string failureCompatible;
    static bool build;
    //    static bool print;
    static unsigned int npos;
    virtual ~Index();
    Index();
    Index(const Index& Other);

    /// Index for product basis (default floor at last basis)
    Index(const std::vector<const BasisAbstract *> Bases, const std::vector<std::string> Names,unsigned int FloorLevel=-1);

    ///@brief construct Index from Stream
    Index(std::istream &Stream, bool Rewind=true /** if true, go to beginning of file, else read from present position */)
        :Index(Stream,{},{},Rewind){}

    /// Exclusively used in the following constructor
    class empty_subtree_exception: std::runtime_error{
    public:
        empty_subtree_exception(): std::runtime_error("Meant to be caught!"){}
    };

    /// new basic contstruct from axes
    Index(const std::vector<Axis> & Ax, const IndexConstraint* Constraint=0,
          std::vector<unsigned int> Pos=std::vector<unsigned int>(),
          std::vector<const Index*> Path=std::vector<const Index*>());

    /// return new Index where grids are replaced with default bases, where possible
    Index* toIndexBasis(std::vector<int> ContractFactor={} /** IF converting to basis, reduce size by factor */) const
    {return toIndexBasis(ContractFactor,std::vector<const Index*>());}


    unsigned int size() const {return _size;}
    const std::string axisName() const;
    const std::string axisSubset() const;
    void setAxisName(std::string Name);
    void setAxisSubset(std::string Subset);
    static const std::string getAxisSubset(std::string AxisName); ///< get subset by name of axisName

    OperatorTree * localOverlap()const; ///< pointer to overlap for present index branch
    OperatorTree * localInvOvr() const; ///< pointer to overlap with all floors inverted
    const OperatorAbstract *overlap() const {return IndexOverlap::get(this);} ///< overlap for present sub-index (if any)
    void setOverlap(const OperatorAbstract *Ovr){IndexOverlap::set(this,Ovr);}  ///< set overlap pointer for this only (not lower hierarchies)
    const Inverse * inverseOverlap() const; ///< inverse overlap for present sub-index (if any)
    void setInverseOverlap(const Inverse *Inv);///< set inverse overlap pointer for this only (not lower hierarchies)

    void assignOverlap(const OperatorAbstract *Ovr); ///< DEPRECATED! assign overlap to this and all possible lower levels

    const BasisAbstract* basis() const;

    void setBasis(const BasisAbstract* Bas);//{basis()=Bas;}

    void testInverseOverlap() const; ///< check error in Ovr*InvOvr and InvOvr*Ovr
    const OperatorAbstract* ovrNonSingular() const;

    void setFloorAuto(std::vector<std::string> & FemAxes); ///<set floor above leaf or at highest FE level
    void setFloor(unsigned int Level); //!< starting from Level, index will be considered contiguous
    void resetFloor(unsigned int Level); //!< move floor to different level
    void cleanFloors(const bool FloorFound=false); //!< remove extra floor levels below first
    unsigned int depthInFloor() const;         //!< depth below floor level

    /// set up the inverse correction, where needed
    void setInverse(Discretization *Disc, const bool Off=false, unsigned int Nskip=0);
    void dvrWeights(std::vector<std::complex<double> > & Weights) const; //!< get the actual dvrWeights as a single vector
    std::vector<const Index*> path() const; ///< return index path from top to present
    unsigned int heightAboveFloor() const; ///< number of hierarchy levels from this to floor
    unsigned int heightAboveBottom() const; ///< number of hierarchy levels from this to bottom

    void nodeCopy(const Index *Node,bool View);
    bool nodeEmpty() const;
    bool nodeEquivalent(const Index *Other) const; ///< two index NODES bases agree
    size_t diagnoseSizeOfNode() const;

    /// two index TREES are equivalent if the have the same structure and the respective nodes are equivalent
    /// OnlyAxes.size()!=0 restricts nodeEquivalent comparison to given axes
    bool treeEquivalent(const Index *Other, const std::vector<std::string> OnlyAxes={}) const;
    bool subEquivalent(std::string Mess="") const; ///< equivalence of all sub-indices (required for tensor-product operators)

    Index* factor(const std::vector<int> Levels, bool Exclude=false) const;

    void leafAdd();

    void purge(unsigned int Height); ///< remove all children of heigh<Height below present, create matching Basis

    unsigned int sizeCompute() const;                                           //!< recompute the size the tree
    inline unsigned int sizeStored() const{return _size;}

    Index* leafAtPos(unsigned int Pos) const; ///< Pos'th leaf of tree
    std::vector<unsigned int> indexAtPos(unsigned int Pos){return leafAtPos(Pos)->index();}
    unsigned int posIndex(const Index *Root=0) const;///< position count relative to Root (unmerged)
    unsigned int posInFloor() const;///< position count to floor node (unmerged)

    unsigned int axisLevel(const Index *ITree) const; ///< level of axis in ITree (=npos if not contained), duplicate names are distinguished by continuity

    unsigned int position(const Index *SubIndex) const;///< position of SubIndex in global count
    unsigned int nSub(const Index *SubIndex) const;///< index of SubIndex in vector<Index*> Idx

    unsigned int globalLength() const; ///< highest global index -lowest global index + 1
    void boundaryIndices(std::vector<unsigned int> &Boundaries) const; ///< create a vector of all global indices shared between elements
    void globalElementBoundary(double ElemBound,std::vector<unsigned int> & GlobalIndex) const; ///< list of global indices matching element boundary
    /// expand global contracted eigenvectors into set of coefficient vectors (AsDual==true interpretes as Dual vectors with overlap multiplied on)
    void unGlobal(UseMatrix & Eigen, bool AsDual, std::vector<Coefficients *> &Evec, std::vector<int> cols) const;
    std::vector<unsigned int> contractedNumbering() const; ///< return map from position in index numbering continuity is contracted
    void contractedNumbering(std::vector<unsigned int> & Numbering, std::vector<unsigned int> & Multi) const; ///< global numbering and multiplicities
    void matrixContract(const UseMatrix &Mat,UseMatrix & GMat) const; ///< compute matrix where continuity indices are contracted
    std::vector<std::complex<double> > dvrWeigContract() const; ///< DVR weights wrt contracted index
    /// wrapper function: get contracted numbering and normalization (for multiplicities)
    void multiplicities(std::vector<unsigned int> &Glob, std::vector<double>& Norm) const;

    std::string hierarchy_no_NONE() const; ///< HACK: do not include final NONE in hierarchy
    std::string hierarchy(unsigned int HybridPath=INT_MAX) const; ///< hierarchy of axis names, if hybrid, choose child(HybridPath)
    std::string coordinates(unsigned int HybridPath=INT_MAX) const; ///< hierarchy without duplicate(=FE) names
    static std::string coordinates(std::string Hierarchy); ///< construct coordinates from hierachy (remove duplicates and sub, spec, etc)

    bool strDataShow(){return axisName()!="NONE";}
    std::string strAxes() const; ///< overview of index structure: (range of) child sizes for each axis name
    std::string strNode(int Level=Tree_defaultKind) const; //!< single node info
    void show(const std::string & text="Index", int Level=0) const; //!< display index
    void axisPlot(std::string File, int Points, double QLow=-DBL_MAX, double QUp=DBL_MAX) const; ///< plot all basis functions on a hierarchy level

    const Index* findAxisStarts(std::string Start) const; ///< first Index where axis name starts with Start
    const Index * firstFloor() const; //!< return floor at the edge of current branch

    Index* axisIndex(const std::string Name) const; //!< index with axis Name

    unsigned int depthOfDuplicate() const; //!< depth() of level with same axisName (=Index::npos - no continuity level)
    unsigned int continuity() const; //!< depth() of continuity level for present index (=Index::npos - no continuity level)
    unsigned int continuity(unsigned int N) const; //! number of levels from present Index to N'th continuity level (=Index::npos - no continuity level)
    unsigned int contDepth() const; //!< number of continuity levels above present

    bool isFem() const;
    bool isAbsorptive() const;
    bool isBottom() const;      ///< bottom: an index node whose next level is/would be all index leafs
    void bottomExpand() const;  ///< expand bottom index by attaching single index leafs (for algorithms that are easier with full Index tree)
    void bottomUnexpand() const;///< undo previous bottom expand
    void bottomExpandAll() const;  ///< expand all bottom indices in tree
    void bottomUnexpandAll() const;///< undo all previous bottom expands
    bool hasFloor() const;
    void setFloor();
    void setKind(const char Kind);
    char indexKind() const {return _indexKind;}
    void unsetFloor();
    bool isRoot() const {return parent()==0;} //!< true root index (criterion: has no parent())
    bool isHybrid() const;

    std::vector<double> grid(std::string Axis) const; ///< returns grid on axis, emtpy if not grid

    const Index * lowerNeighbor() const;  //!< pointer to lower neighbor of Index, =0 if none
    const Index * upperNeighbor() const;  //!< pointer to upper neighbor of Index, =0 if none
    const Index * lowerNeighbor(unsigned int D) const;  //!< pointer to lower neighbor of leaf in D-direction, error if not leaf
    const Index * upperNeighbor(unsigned int D) const;  //!< pointer to upper neighbor of leaf in D-direction, error if not leaf

    void localOverlapAndInverse(OperatorTree *Ovr, OperatorTree *Inv);                   //!<set up a default overlap for Index
    void getS0inv(UseMatrix & S0inv,unsigned int D) const; //! tensor factor of inverse overlap in D-direction (error if not tensor factor)

    void writeStructure(std::ofstream & stream) const;            ///< write index structure to file
    bool compatibleFile(std::ifstream & stream,int code=0) const; ///< true if structure on file is compatible with index

    ///@brief write full definition of Index to present position of Stream
    void write(std::ostream & Stream, bool Enter=true /** internal use, do not set */ ) const;
    bool isOverlapDiagonal() const;

    /// Convenience wrapper around BasisAbstract::physical
    double physical() const;

    /// establishes unambiguous match with Other index (requires dynamic_cast to derived classes)
    virtual const Index* from() const;

};

#endif
