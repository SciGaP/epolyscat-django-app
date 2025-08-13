// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORSINGLE_H
#define OPERATORSINGLE_H

#include "tools.h"
#include <deque>
#include "useMatrix.h"
#include "labelled.h"

class Discretization;
class Index;
//class IndexFloor;
class Coefficients;
class CoefficientsFloor;

/** \ingroup OperatorFloors */
/// (OBSOLETE)
class OperatorSingle : public Labelled<OperatorSingle> {
//    friend class Operator;
    friend class Coefficients;
public:

    static std::deque<UseMatrix> matsTable; /// table to contain all mats pointers of operator single
    static std::map<std::string,UseMatrix> matsTableNew; /// table to contain all mats pointers of operator single
    static const UseMatrix *matsAdd(UseMatrix &mat,std::string Hash="noHash"); //the argument 'mat' can (and should) be deleted after calling this function

    void static fuse(std::vector<OperatorSingle *> &Ops);

    std::string name;
    std::string definition;
//CHANGE    std::complex<double> * factor; // SHOULD GO PRIVATE apply matrices and multiply by factor
    const Index *iIndex,*jIndex; ///< left/right floor indices
//    OperatorSingle(const Index *iIndex = 0,const Index *jIndex = 0);

    virtual ~OperatorSingle();

    OperatorSingle(const std::string Name="noName", const std::string Def="unDef", const Index *IIndex=0, const Index *JIndex=0);
    std::vector<const UseMatrix*> mats;

    // === functions =========================================
    virtual void axpy(std::complex<double> Alfa, CoefficientsFloor & X,std::complex<double> Beta, CoefficientsFloor & Y, bool transpose=false) const=0; //Y = Operator*X+const;
    void axpy(CoefficientsFloor & X, CoefficientsFloor & Y, bool transpose=false) const {axpy(1.,X,1.,Y,transpose);} //Y = Operator*X+const;

//    virtual void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
//                      const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const;

    virtual void apply(std::complex<double>* InOut) const; //new version for floors Y = Operator*X+const;
    virtual void apply(std::vector<std::complex<double> > & InOut) const; //Y = Operator*X+const;
    virtual bool isZero(double Eps=0.) const=0;
    virtual void inverse() = 0; ///< replace operator with its inverse
    virtual void matrix(UseMatrix & mat) const; ///< set up matrix

    virtual std::string str() const {return name;}
    std::string strStructure() const;
    virtual std::string strDataStructure() const{return "";} //!< supplement with specific data information

    void show(std::string Mess){if(Mess!="")std::cout<<Mess+": "; std::cout<<strStructure()<<std::endl;}
    double norm() const;

//    virtual void packInfo(std::string &Info);
//    virtual void unpackInfo(std::string &Info);

protected:
    virtual double setNorm(); //!< infty-Norm
    double oNorm;
    double oNonzeros;
    virtual bool absorb(OperatorSingle * & Other){return false;}
    void showStructure(const std::string Text="") const {std::cout<<Text<<": "<<strStructure()<<std::endl;} //!< print structure information
    const Discretization *dataDisc(const Discretization *iDisc, const Discretization *jDisc); // substitutes askjdsisc. // points to where to get data
    std::complex<double> * timeDepFac;
};


#endif // OPERATORSINGLE_H
