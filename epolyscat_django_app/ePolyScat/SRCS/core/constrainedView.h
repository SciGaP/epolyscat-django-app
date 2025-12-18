// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <string>
#include <vector>
#include <memory>

#include "tree.h"
#include "operatorTree.h"
#include "str.h"
#include "operatorFloor.h"
#include "projectSubspace.h"

class OperatorFloor;
class OperatorTree;

///\brief constrained view of on an operator
class ConstrainedView: public OperatorTree
{
protected:
    const OperatorAbstract * _op;
public:
    static const OperatorAbstract* factory(OperatorTree *Op, std::string Constraint="");

    virtual void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const=0;

    const OperatorAbstract * derivedFrom() const {return _op;}
    virtual bool isBlockDiagonal() const;

    std::string strNode() const;
};

class ConstrainedZero: public ConstrainedView
{
    friend class OperatorFloor;
    void addZeros(std::string Axis,const std::vector<double> & Zeros,const Index* Idx,std::vector<unsigned int> & IZeros);
    class Floor: public OperatorFloor{
        friend class ConstrainedZero;
        OperatorFloor *_parentFloor;
        std::vector<unsigned int> _rightZeros,_leftZeros;
    public:
        Floor(OperatorFloor* Parent,std::vector<unsigned int>LeftZeros,std::vector<unsigned int>&RightZeros);
        Floor(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf);

        void axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX, const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const;
        void pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const;
    };


public:
    ConstrainedZero(OperatorTree* Op,std::vector<std::string> Axes, const std::vector<std::vector<double> > & Zeros);
    bool isView() const { return _view; }
    bool isDummy() const;
    const OperatorFloor* floor() const override;
    OperatorFloor*& floor() override;
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
        OperatorTree::apply(A,Vec,B,Y);
    }
};

class ConstrainedZeroDerivative: public ConstrainedView{
    void getTrans(std::string Axis, const std::vector<double> & Zeros, const Index* Idx,
                  std::vector<std::vector<std::complex<double> > > &UTrans);
    class Floor: public OperatorFloor{
        OperatorFloor *_parentFloor;
        std::vector<std::vector<std::complex<double> > > _uRight,_uLeft;
    public:
        Floor(OperatorFloor* Parent,std::vector<std::vector<std::complex<double> > > & ULeft,
              std::vector<std::vector<std::complex<double> > > & URight);
        void axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX, const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const;
        void pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const{DEVABORT("to be implemented");}
        void project(std::complex<double>*Vec,int Size, const std::vector<std::vector<std::complex<double> > > & U) const;
    };
public:
    ConstrainedZeroDerivative(OperatorTree* Op,std::vector<std::string> Axes, const std::vector<std::vector<double> > & Zeros);
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
    {OperatorTree::apply(A,Vec,B,Y);}
};
#endif // CONSTRAINT_H
