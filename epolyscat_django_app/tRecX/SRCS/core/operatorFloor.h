// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef OPERATORFLOOR_H
#define OPERATORFLOOR_H

#include <map>
#include <string>
#include <vector>
#include <complex>
#include <memory>
#include "abort.h"

#include "qtEigenDense.h"
#include "coefficients.h"

class UseMatrix;
class Index;
class OperatorTensor;
class OperatorDefinition;

/** @defgroup OperatorFloors Lowest level
 *  \ingroup Operators
 *  \brief Fast and data-contigous: full matrix, diagonal, sparse, tensor products, specialized forms, etc.
 *  @{
*/

/// \brief lowest level operator
class OperatorFloor
{
    //HACK
    friend class HaccInverse;
    friend class OperatorFloorSpecial;
    friend class Parallel;

    static std::string tensorType(std::string Kinds);
    static std::map<std::string,int> _packCode;
protected:
    static std::map <std::string,std::vector<std::complex<double> > > complexData;
    static std::map <std::string,std::vector<double> > realData;
    static std::vector<std::complex<double> > tempComplex;
    static double bandedRatio; // ratio of non-zeros when banded is used
    static std::string hashString(unsigned int Rows, unsigned int Cols);

    static void scale(std::complex<double> Beta, std::vector<std::complex<double> > & Y);
    static void scale(std::complex<double> Beta,std::complex<double> *Y, unsigned int SizY);
    static UseMatrix UseMatrixTensor(std::vector<const UseMatrix*> Dat);

    std::shared_ptr<std::vector<std::complex<double> > > addComplex(const std::string & Hash, const std::vector<std::complex<double> > &Dat);
    std::vector<double > * addReal(const std::string & Hash, const std::vector<double> &Dat);

    std::string _kind;
    std::shared_ptr<std::vector<std::complex<double>>> iWeights;
    unsigned int _rows,_cols;
    double oNorm;
//    std::vector<std::complex<double> > * dat;
    std::shared_ptr<std::vector<std::complex<double> >>  dat;

    /// generic norm calculator - can be very slow, overwrite by specific wherever possible
    virtual void setNorm() const;

    /// return basic information about the OperatorFloor
    /// - size of specific buffer (to be filled by pack())
    /// - packCode                (to be inserted by pack())
    /// - rows
    /// - cols
    /// - iWeight->size()
    /// - Buf contains *iWeights
    void packBasic(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const;
    void unpackBasic(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf);

    double _cost;
    std::complex<double> * timeDepFac;

    virtual void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
                      const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const=0;

    virtual void axpyTranspose(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
                               const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const
    {ABORT("not implemented for "+kind());}
public:
    static int diagnoseDataSize(); ///< current total floor data size

    static double UNDEFINED;
    virtual ~OperatorFloor();
    OperatorFloor(std::string Kind):_kind(Kind),iWeights(0),_rows(0),_cols(0),oNorm(UNDEFINED),dat(0),_cost(UNDEFINED),timeDepFac(0){}
    OperatorFloor(unsigned int Rows,unsigned int Cols,std::string Kind):_kind(Kind),iWeights(0),_rows(Rows),_cols(Cols),oNorm(UNDEFINED),_cost(UNDEFINED),timeDepFac(0){}

    /// new OperatorFloor from PMats, Hash is used for data registry
    static OperatorFloor* factory(const std::vector<const UseMatrix *> &PMats, std::string Hash);

    /// new OperatorFloor from PMats, Hash is used for data registry (Eigen version)
    static OperatorFloor* factory(const std::vector<const Eigen::MatrixXcd*> &PMats, std::string Hash);

    /// new OperatorFloor from OperatorDefintion
    static OperatorFloor* factory(const std::string& TermOper,
                                  const Index * IIndex,const Index *JIndex, std::complex<double> Multiplier);

    /// new OperatorFloor from OperatorTensor
    static OperatorFloor* factory(const OperatorTensor &OT);

    /// new OperatorFloor from OperatorFloor
    static OperatorFloor* copyFactory(const OperatorFloor *Other);

    /// new banded "inverse" OperatorFloor from Idx; needed for CoulX
    static OperatorFloor* factoryInverse(const Index* Idx, unsigned int SubD, unsigned int SuperD, bool BandOvr);

    /// new OperatorFloor from Info and Buf
    static OperatorFloor * unpackFactory(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf);

    /// factory for definition of special [[...]]-operators
    static OperatorFloor* specialFactory(const std::string Name, const std::string Def,
                                         const Index *IIndex, const Index *JIndex,std::complex<double> Multiplier);

    /// replace Floor with new floor constructed from Info and Buf (empty - construct dummy)
    static void replace(OperatorFloor *&Floor,
                        const std::vector<int> & Info=std::vector<int>(),
                        const std::vector<std::complex<double> >&Buf=std::vector<std::complex<double> >());

    static bool absorb(OperatorFloor *&Into, OperatorFloor *&Other, std::string Def);

    /// apply floor including time-dependent factors
    virtual void apply(std::complex<double> Alfa, const std::complex<double>* X, unsigned int SizX,
                       const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const;

    virtual void matrix(UseMatrix & mat) const; ///< set up matrix
    virtual Eigen::MatrixXcd matrixFactor(int D) const; ///< D'th matrix factor (=matrix for non-tensor)
    virtual Eigen::MatrixXcd matrix() const; ///< return floor matrix

    virtual double norm() const{if(oNorm==UNDEFINED)setNorm();return oNorm;}
    virtual bool isAbsorbable() const;


    std::string str(int Digits=0) const;
    virtual std::string strInfo() const; ///< for neat printing
    void write(std::ofstream &File) const;
    static OperatorFloor* readFactory(std::ifstream &File);

    virtual double applicationCost(bool Bcast=true) const;

    virtual bool isZero(double Eps=0.) const {return oNorm<=Eps;}
    virtual bool isDiagonal() const;
    void applyLeftOverlap(std::complex<double> * Data);
    void addWeights(const Index *IIndex,const std::vector<std::complex<double> > & IWeights);

    static int packCode(std::string Kind);
    virtual void pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const=0;

    virtual void uniqueData(){ABORT("not implemented for "+kind());} ///< make sure OperatorFloor does not share data storage
    std::complex<double> * factor() const {return const_cast<OperatorFloor*>(this)->timeDepFac;}
    void setFactor(std::complex<double>* Factor) {timeDepFac=Factor;}
    std::string kind() const {return _kind;}

    static std::string failAbsorb; // for debugging
    unsigned int rows() const {return _rows;}
    unsigned int cols() const {return _cols;}
    
    virtual long applyCount() const{ ABORT("applyCount not implemented"); }

    /// broad-cast non-zero Ofloor to all other threads
    static OperatorFloor* bCast(const OperatorFloor* Ofloor);

    /// forced reset of Norm
    void forceNorm(double Norm=UNDEFINED) const {const_cast<OperatorFloor*>(this)->oNorm=Norm;}
    /// forece reset of application cost
    void forceCost(double Cost=UNDEFINED) const {const_cast<OperatorFloor*>(this)->_cost=Cost;}

};

/// placeholder dummy floor
class OperatorDUM: public OperatorFloor{
public:
    OperatorDUM(double Norm=UNDEFINED):OperatorFloor(0,0,"DUM"){oNorm=Norm;} // can be constructed with a dummy norm
    void pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const {packBasic(Info,Buf);}
    void axpy(const std::complex<double> & Alfa, const std::complex<double>*X, unsigned int SizX,
              const std::complex<double> & Beta,       std::complex<double>*Y, unsigned int SizY) const;
    double applicationCost(bool Bcast=true) const {if(_cost==UNDEFINED)ABORT("OperatorDUM has not been assigned an application cost");return _cost;}
};
/// 1x1 single number floor
class Operator1X1: public OperatorFloor{
    std::complex<double> _c;
public:
    Operator1X1(std::complex<double> C):OperatorFloor(1,1,"1X1"),_c(C){oNorm=std::abs(C);} // can be constructed with a dummy norm
    void pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const {packBasic(Info,Buf);}
    void axpy(const std::complex<double> & Alfa, const std::complex<double>*X, unsigned int SizX,
              const std::complex<double> & Beta,       std::complex<double>*Y, unsigned int SizY) const {*Y=Alfa*_c*(*X)+Beta*(*Y);};
};
/// placeholder Zero floor
class OperatorZero: public OperatorFloor{
public:
    OperatorZero(int Rows, int Cols):OperatorFloor(Rows,Cols,"Zero"){oNorm=0.;}
    void pack(std::vector<int> &Info, std::vector<std::complex<double> >&Buf) const {packBasic(Info,Buf);}
    void axpy(const std::complex<double> & Alfa, const std::complex<double>*X, unsigned int SizX,
              const std::complex<double> & Beta,       std::complex<double>*Y, unsigned int SizY) const;
    double applicationCost(bool Bcast=true) const {return 0.;}
};

/** @} */ // end group OperatorFloors

#endif // OPERATORAPPLY_H
