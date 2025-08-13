// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INVERSEFEM_H
#define INVERSEFEM_H

#include "inverse.h"
#include "operatorTree.h"
#include "operatorFloor.h"
#include "mpiWrapper.h"
#include "parallelContinuity.h"

class Coefficients;
class OperatorTree;
class ParallelContinuity;

class InverseFEM : public Inverse
{
    static std::map<std::string,ParallelContinuity*> _contTab;

    class Floor: public OperatorFloor{
        const Index* iIndex,*jIndex;
        int _marginDepth; // depth of margin index (relative to IdxRoot)
        bool _upperMargin;
        std::vector<std::complex<double> > _s0InvM; // the correction vector
        void axpyRecursive(const Index* XIdx,const Index* YIdx, const std::complex<double> & Alfa, const std::complex<double> *&C, const std::complex<double>* &X, std::complex<double> *&Y) const;
   public:
        Floor(Coefficients* S0InvM, Coefficients* RtS0invR, std::complex<double> Multi);
        void axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX, const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const;
        void pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const{DEVABORT("not implemented");}
    };

    class CorrectionMap: public OperatorTree{
    public:
        CorrectionMap(Coefficients* S0InvM, Coefficients* Marg, std::complex<double> Multi);
    };
    OperatorTree * _correctionMap;

    int _nSplit;
    Coefficients* _margin;
    ParallelContinuity * _contMargin; // margin pointers for _margin;

    void constructCorrection(Index* Idx);

public:
    ~InverseFEM(){delete _correctionMap;delete _margin; delete _contMargin;}
    InverseFEM(Index *Idx, int Begin=-1, int End=-1);
    void parallelSetup() const {if(MPIwrapper::Size()>1)ABORT("cannot be run in parallel");}
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    void apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{DEVABORT("not implemented");}
    void applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;

    void apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{if(MPIwrapper::Size()>1)DEVABORT("cannot be run in parallel");}
    void applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const{if(MPIwrapper::Size()>1)DEVABORT("cannot be run in parallel");}

};


#endif // INVERSEFEM_H
