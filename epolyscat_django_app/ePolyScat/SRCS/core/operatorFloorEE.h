// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORFLOOREE_H
#define OPERATORFLOOREE_H

#include "operatorFloor.h"
#include "basicDisc.h"
#include "readInput.h"
#include "useMatrix.h"

using namespace std;

/** \ingroup OperatorFloors */
/// \brief two-particle Coulomb repulsion, using multipole expansion
class OperatorFloorEE : public OperatorFloor
{
    std::vector<int > Idx,Jdx;
    void construct_index(const Index* I, std::vector<int> &Idx); // order n1,n2,m1,m2,l1,l2
     std::vector<std::complex<double> > qWeights;

    // Helpers //////////////
    class BasicDiscOnlyRadial: public BasicDisc{
    public:
        BasicDiscOnlyRadial(ReadInput &In);
        int lmax;
    };

    class Helpers_e_e{
        static void lobatto_quadratureNew(const BasisAbstract*, Eigen::VectorXd&, Eigen::VectorXd&);
        static void analyseDVR(const BasisAbstract* B, UseMatrix& pts, UseMatrix& wgs, vector<bool>& wgs_inc);
    public:
        static std::vector<std::vector<std::vector<UseMatrix > > > TwoElecOverlap;// TwoElecOverlap[n1][n2][lambda]
        static std::vector<std::vector<std::vector<double> > > TwoElecOverlapMaxValue;
        static void initialize(const Index *Root, int& lambda_upper_limit);
        static void absorbInverse(const Index *Root);
    
    };
    // /////////////////////
    static double epsGaunt;
    static double epsTwoElec;

    static int lambda_upper_limit;
protected:
    void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
              const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const;

    /// OBSOLESCENT
    void axpy(std::complex<double> Alfa, const std::vector<std::complex<double> > & X,
                      std::complex<double> Beta, std::vector<std::complex<double> > & Y) const;
public:
    static int read(ReadInput& Inp);
    static void constrain(UseMatrix& Mult, const Index* IIndex, const Index* JIndex);

    OperatorFloorEE(const std::string Name, const std::string Def, const Index *IIndex, const Index *JIndex);
    OperatorFloorEE(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf);
    ~OperatorFloorEE();
    bool isZero(double Eps=0.) const;

    void pack(std::vector<int> & Info,std::vector<std::complex<double> > &Buf) const;
};

#endif // OPERATORFLOOREE_H
