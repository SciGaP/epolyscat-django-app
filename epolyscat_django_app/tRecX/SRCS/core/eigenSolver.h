// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <vector>
#include <string>
#include "tree.h"
#include "eigenSolverAbstract.h"

class OperatorAbstract;
class Index;
class Resolvent;

#include "operatorTree.h"

#include "arpack.h"
#include "coefficients.h"
#include "coefficientsSparse.h"

/// \ingroup Linalg
/// \brief Eigen solver for OperatorAbstract - detects and exploits block-structure
///
/// optionally computes dual of right-eigenvectors
///
/// the projector onto the subspace of energies is
///
/// P = sum[i] rightEigenvector[i] dualVector[i].transpose()
///
class EigenSolver: public EigenSolverAbstract, public Tree<EigenSolver>
{   
    std::string _method;
    bool _fullVectors;

    void expandInto(std::vector<Coefficients *> &CVec, std::vector<std::complex<double> > &Val,
                                 const std::complex<double> *PVal,  const Eigen::SparseMatrix<std::complex<double>> & Vec, double Emin, double Emax,
                                 bool ExcludeRange);
    void _compute();
    void _computeBlock(const OperatorTree *OpTree, const OperatorAbstract *Ov);


    ///\brief fix for unique phases and normalization (if ambiguous)
    ///
    /// general non-symmetric: sum_i |(R_m)_i|^2=1 for all m
    /// hermitian: coefficient at maximal component >0

    void collect(); // establish eigenvalues and eigenvectors on present level

public:
    virtual ~EigenSolver();

    /// defaults to all eigenvalues and right eigenvectors
    EigenSolver();

    /// compute eigenvalues in range
    EigenSolver(double Emin /** lowest eigenvalue real part to include */,
                double Emax /** highest eigenvalue real part to include */,
                bool RightVectors=true /** compute right eigenvectors */,
                bool DualVectors=false /** compute dual vectors */,
                bool ExcludeRange=false /** compute eigenvalues OUTSIDE [Emin,Emax] */,
                std::string Method="auto" /** Lapack, Arpack */
            )
        :EigenSolver(Emin,Emax,INT_MAX,RightVectors,DualVectors,ExcludeRange,Method){}

    EigenSolver(double Emin,double Emax,int Nmax,bool RightVectors, bool DualVectors, bool ExcludeRange, std::string Method);


    std::vector<std::complex<double> > eigenvalues(); ///< collect eigenvalues from all blocks

    // Does NOT retain ownership of these
    EigenSolver & fullVectors();//{_fullVectors=true;return *this;}
    std::vector<Coefficients*> rightVectors(); ///< right hand eigenvectors right_j
    std::vector<Coefficients*> dualVectors(); ///< duals of eigenvectors:  dual_i:=(left_i)^* S,  = dual_i.inner(right_j)=delta_ij



private:

    class Arp:public Arpack{
        static std::complex<double> shift; // subtract shift*S from Hamiltonian to move desired eigenvalues to below 0
    protected:
        // internal data
        const OperatorAbstract* o;
        std::vector<std::complex<double> *> pX,pY;
        Coefficients x,y,tmp;
        void apply(const std::complex<double> *X, std::complex<double> *Y);
    public:
        std::complex<double> & ref_shift(){return shift;}
        Arp(const OperatorAbstract* O);
        virtual void eigen(std::vector<std::complex<double> > &Eval, std::vector<Coefficients* > &Rvec,
                           unsigned int Nvec, const std::string &Which, bool Restart);
    };

    class ArpShiftInverse: public Arp {
        std::unique_ptr<Resolvent> _resolv;
        void apply(const std::complex<double> *X, std::complex<double> *Y);
    public:
        ArpShiftInverse(const OperatorAbstract* O, const std::complex<double> Eguess);
        void eigen(std::vector<std::complex<double> > &Eval, std::vector<Coefficients* > &Rvec,
                   unsigned int Nvec, const std::string &Which, bool Restart);
    };

};

#endif // EIGENSOLVER_H
