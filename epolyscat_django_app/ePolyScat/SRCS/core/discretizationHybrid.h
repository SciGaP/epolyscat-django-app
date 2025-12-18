// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONHYBRID_H
#define DISCRETIZATIONHYBRID_H

#include "discretization.h"
#include "operatorTree.h"
#include "inverse.h"
#include "projectSubspace.h"

class BasicDisc;
class OperatorVectors;

/** \ingroup Discretizations */

///@brief Direct sum of several BasicDisc's
///
/// usage examples:
///
/// 210OffCenter.inp,220HybridSubspace.inp
class DiscretizationHybrid : public Discretization
{
//    friend class Operator;
    friend class OperatorTree;
    friend class DiscretizationSurface;
    friend class DiscretizationSurfaceHybrid;
    friend class OperatorMapChannelsSurface;
    friend class Evaluator;
    friend class BasisNdim;

    std::vector<std::string> compName;
    std::vector<const BasicDisc*> comp;

    // a new class should enclude this, Hybrid should be abstracted
    std::unique_ptr<ProjectSubspace> _projSub;

    class Overlap:public OperatorTree
    {
        std::vector<std::vector<const OperatorTree*> > _block;
    public:
        Overlap(const DiscretizationHybrid  * H, std::string OffDiag);
        void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
        const OperatorTree * block(int K, int L) const {return _block[K][L];}
    };

    class Sinverse:public Inverse
    {
        OperatorVectors *sInvBA;   // Sb^-1 C
        const OperatorAbstract *sAB;   // C^H
        const OperatorAbstract *sAinv,*sBinv;
        Eigen::MatrixXcd zInv; // Z^-1
        Coefficients *aVec,*bVec;
        Sinverse(const Sinverse & Other); //do not allow copy        
    public:
        ~Sinverse();
        Sinverse(DiscretizationHybrid* H, const Overlap* Ovr);
        void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
        void apply(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const;

        void applyCorrection(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
        {ABORT("to be implemented");}
        void apply0(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const
        {ABORT("to be implemented");}
        void applyCorrection(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const
        {ABORT("to be implemented");}
        void apply0(std::complex<double> A, const CoefficientsLocal &Vec, std::complex<double> B, CoefficientsLocal &Y) const
        {ABORT("to be implemented");}
        void parallelSetup() const;
    };

public:
    static bool isHybrid(ReadInput & Inp, int Line=1);
    /// compose hybrid operator from non-hybrid blocks
    DiscretizationHybrid(ReadInput &In); ///< Construct using input file
    DiscretizationHybrid(const Index* Idx); ///< Temporary: construct from existing Index
    const Discretization* haCC() const;
};

#endif // DISCRETIZATIONHYBRID_H
