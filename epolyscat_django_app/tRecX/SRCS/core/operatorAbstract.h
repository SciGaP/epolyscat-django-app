// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef OPERATORABSTRACT_H
#define OPERATORABSTRACT_H
#include "abort.h"
#include <vector>
#include <complex>
#include <cfloat>
#include <memory>
#include "linSpaceMap.h"
#include "labelled.h"

#include "qtEigenDense.h"
#include "qtEigenSparse.h"
#include "tree.h"

#include "operatorDefinition.h"

class Coefficients;
class Index;
class UseMatrix;
class Plot;
class ReadInput;

/** @defgroup Operators
 * \brief definitions, structures, evaluation
*/


/** @defgroup Structures
 * \ingroup Operators
 * \brief maps between Coefficients: recursive, spectral projections, diagonal operators, etc.
 *  @{
*/

/// Abstract base class for any map between recursively Index'd spaced
class OperatorAbstract: public LinSpaceMap<Coefficients>, public Labelled<OperatorAbstract>{
    mutable Coefficients * _tempRHS;
    mutable Coefficients * _tempLHS;

protected:
    double _time;
    /// add contract overall matrix and add into GMat (optionally times Factor)
    virtual void matrixContract(const UseMatrix &Mat, UseMatrix &GMat, std::complex<double> Factor=1.,
                                std::vector<int>ISort=std::vector<int>(0),std::vector<int>JSort=std::vector<int>(0)) const;
    OperatorDefinition definition;
    virtual bool isSymmetric(std:: string  Kind, double Eps=1.e-12) const;

public:
    static bool useOperatorFloor;
    static bool fuseOp;
    static bool useTensor;
    static bool flat;
    static std::string eigenMethod;
    static std::complex<double> Arp_shift;
    static void readControls(ReadInput &Inp);

    std::string name;
    const Index *iIndex, *jIndex;     ///< left and right indices of operator

    OperatorDefinition def() const {return definition;}

    virtual ~OperatorAbstract();

    /// a very simple factory
    static std::shared_ptr<OperatorAbstract> factory(std::string Name,std::string Definition,const Index* Idx, const Index* Jdx);

    OperatorAbstract():_tempRHS(nullptr),_tempLHS(nullptr),name("undefined"),iIndex(nullptr),jIndex(nullptr){}
    OperatorAbstract(std::string Name,const Index* IIndex,const Index* JIndex)
        :_tempRHS(nullptr),_tempLHS(nullptr),_time(0.),name(Name),iIndex(IIndex),jIndex(JIndex){}

    void setDefinition(std::string Def){definition=OperatorDefinition(Def);}

    const Index* idx(bool True) const{return True?iIndex:jIndex;} ///< left (True=true) or right index
    Index* idx(bool True); ///< left (True=true) or right index
    const Index* idx() const{return iIndex;}
    const Index* jdx() const{return jIndex;}
    Index* idx() {return const_cast<Index*>(iIndex);}
    Index* jdx() {return const_cast<Index*>(jIndex);}


    /// required for abstract base class LineSpaceMap
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const=0;
    //    virtual void update(double Time, const Coefficients* CurrentVec=nullptr); ///< required for LineSpaceMap
    //    const OperatorAbstract& time(double Time){update(Time);return *this;} ///< update time and return operator
    const Coefficients & lhsVector() const {return *tempLHS();}///< required for LineSpaceMap
    const Coefficients & rhsVector() const {return *tempRHS();}///< required for LineSpaceMap

    virtual void axpy(std::complex<double> Alfa,const Coefficients& X,std::complex<double> Beta, Coefficients& Y, const double Time);

    virtual void axpy(std::complex<double> A, const Coefficients &X, std::complex<double> B, Coefficients &Y){apply(A,X,B,Y);}
    virtual void axpy(const Coefficients &X,Coefficients &Y){apply(1.,X,1.,Y);}

    Coefficients* tempLHS() const; ///< pointer to storage of lhs vector
    Coefficients* tempRHS() const; ///< pointer to storage of rhs vector

    /// print minimal info
    virtual std::string str(int Level=Tree_defaultKind) const;// {return name+" (no structure information defined)";}

    /// convert to full matrix, column-wise storage in Mat (slow)
    virtual Eigen::MatrixXcd matrix() const; ///< convert to full matrix

    /// sparse matrix for operator, near-zeros by Eps are omitted (see BlockView for exact interpretation of Eps)
    virtual Eigen::SparseMatrix<std::complex<double>> matrixSparse(bool Contract=false /** true: contract w continuity conditions */, double Eps=1.e-12) const;

    virtual void matrix(std::vector<std::complex<double> > & Mat) const;
    virtual void matrix(UseMatrix & Mat) const;///< convert to full matrix, column-wise storage in Mat (slow)

    virtual double norm() const;

    UseMatrix matrix(std::vector<Coefficients*> Left, std::vector<Coefficients*> Right) const;

    /// convert to sub matrix, column-wise storage in Mat (slow)
    virtual void subMatrix(std::vector<std::complex<double> > & Mat,const Index* ISub, const Index* JSub) const;

    /// add operator to matrix, create new if empty matrix
    virtual void matrixAdd(std::complex<double> factor, UseMatrix & Mat) const;

    /// get matrix boundary conditions are imposed on coeffecicitens
    ///
    /// dimension is reduced compared to uncontracted
    virtual Eigen::MatrixXcd matrixContracted() const;

    /// (pseudo-)orthonormalos wrt to Op: 1=<C|Op C> (<C^*|Op C> for pseudoScalar)
    void orthoNormalize(std::vector<Coefficients *> &C, bool pseudo=false) const;

    /// (pseudo-)orthonormalize near-degenerate subspaces wrt. Op
    void orthoNormalizeDegenerate(std::vector<std::complex<double> > Eval, std::vector<Coefficients *> Evec, bool Pseudo) const;

    /// <wf1|Op|wf2>, if pseudoScalar <wf1*|Op|wf2>
    virtual std::complex<double> matrixElement(const Coefficients &Ci, const Coefficients &Cj, bool pseudoScalar = false) const;

    std::complex<double> matrixElementUnscaled(const Coefficients &Ci, const Coefficients &Cj) const;///< <wf1*|Op|wf2>
    
    double applicationCost() const;
    virtual long applyCount() const{ ABORT("applyCount not implemented"); }

    virtual bool isBlockDiagonal() const;///< only Tree-type operators can be block-diagonal

    virtual bool isComplexSymmetric(double Eps=1e-12) const; ///< true if M(i,j)=M(j,i) for all i,j (Not necessarily hermitian!)
    virtual bool isSelfAdjoint(double Eps=1e-12) const; ///< true if conj(M(i,j))=M(j,i) for all i,j

    virtual OperatorDefinition getDefininition() const {return definition;}

    bool isHuge(std::string Message="") const;
    virtual bool isIdentity(double Eps=1.e-12, bool Stochastic=true) const;
    virtual bool isZero(double Eps=0.,bool Stochastic=true) const;
    virtual bool isDiagonal() const;

    std::ofstream* write(std::ofstream *File=0); ///< write name, iIndex, jIndex, create default File if needed and return
    std::string plotActionOnCoefficients(const Coefficients * C=0, std::shared_ptr<Plot> Plt=0) const;
    OperatorAbstract(std::ifstream &File, std::string Name, const Index *IIndex, const Index *JIndex); ///< create with name, iIndex, jIndex from File
};
/** @} */ // end group Operators


#endif // OPERATORABSTRACT_H
