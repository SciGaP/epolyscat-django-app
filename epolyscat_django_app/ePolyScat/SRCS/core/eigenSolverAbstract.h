// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef EIGEN_SOLVER_ABSTRACT_H
#define EIGEN_SOLVER_ABSTRACT_H

#include <vector>
#include <complex>
#include <string>

#include "operatorAbstract.h"
#include "index.h"
#include "abort.h"
#include "str.h"


class Coefficients;
class ReadInput;

/// \ingroup Linalg
/// \brief Abstract base class for eigensolver on OperatorAbstract
///
/// Solves the generalized linear eigenproblem (possibly on a subspace)
/// <br>Op vec = S vec eval
/// <br> S is invertible on a subspace
/// <br> Op is an OperatorAbstract and S is provided through the Operator's Index
///
/// the class holds configuration info, which can be altered by various setters <br>
/// construct / call compute() / retrieve results such as eigenvalues

class EigenSolverAbstract {
    std::unique_ptr<OperatorTree> _locOvr;
    void _postProcess();
protected:
    bool _serial; ///< do not use MPI parallelism
    const OperatorAbstract* _op;
    const OperatorAbstract * _ovr;

    std::vector<std::complex<double> > _eigenvalues;
    std::vector<Coefficients* > _rightVectors;
    std::vector<Coefficients* > _leftVectors;
    std::vector<Coefficients* > _dualVectors;

    bool _computeLeftVectors;
    bool _computeRightVectors;

    std::string _sort;  ///< how to sort: SmallReal, SmallAbs, LargeReal, LargeAbs etc. see tools::sortKey for all options

    mutable double _epsVerified; ///< keeps track to which level current solutions have been verified (avoid repeated verification)
    virtual void _compute()=0; ///< solve the raw eigenproblem (w/o sorting, orthonormalization, etc.)

    /// determine a reasonable eigensolver method
    std::string autoMethod(const std::string Method);

    std::string _select;
    /// generate a select string from old-style input
    void makeSelect(int Nev, double Emin, double Emax, bool ExcludeRange);
    int numberEigenv() const; ///< desired number of eigenvalues
    double eMax() const;      ///< desired maximal eigenvalue (real part)
    double eMin() const;      ///< desired minimal eigenvalue (real part)
    bool excludeRange() const;///< true: ask for eigenvalues outside intervals
    std::complex<double> guessEigenvalue() const;
public:

    static void readControls(ReadInput & Inp);
    static std::string method,defaultMethod;
    static void setMethod(std::string Method); ///< set method (Lapack, Arpack...) for eigensolver
    static void resetMethod(){method=defaultMethod;}///< reset method to defaultMethod
    /// Yet to be implemented in a useful way
    static EigenSolverAbstract* factory(OperatorAbstract* _op);

    virtual ~EigenSolverAbstract();
    EigenSolverAbstract(const OperatorAbstract* Op=0);
    void withSelect(std::string Select); ///< update the eigenvalue selection string with Select
    std::string selection() const {return _select;}

    void computeLeftVectors(bool Compute){ _computeLeftVectors=Compute; }
    void computeRightVectors(bool Compute){ _computeRightVectors=Compute; }


    virtual std::vector<std::complex<double> > eigenvalues();
    virtual std::vector<Coefficients* > leftVectors(); ///< left hand eigenvectors: (left_i)^* A = e_i (left_i)^* S
    virtual std::vector<Coefficients* > rightVectors(); ///< right hand eigenvectors Op * right_j = S * right_j eval_j
    virtual std::vector<Coefficients* > dualVectors(); ///< dual vectors (dual_i).innerProduct(right_j) = delta_ij

    bool isComplexSymmetric() const; ///< Complex symmetric eigenproblem A=A^T and S=S^T
    bool isSelfAdoint() const;       ///< Selfadjoint eigenproblem A=A^* and S=S^*

    /// generate duals for symmetric or selfadjoint problems
    std::vector<Coefficients* > symmetricDuals();

    ///\brief orthognalize eigenvectors (and matching duals) in degenerated subspaces
    void orthonormalizeDegenerate(double Eps /** if |E1-E2|<Eps E1,E2 are degenerated */,
                                  const std::vector<std::complex<double> > & Eval,const std::vector<Coefficients*> & Dual,
                                  const std::vector<Coefficients*> & Rvec);

    ///\brief Schmidt-style orthonormalize right eigenvectors and duals
    ///
    /// Dual[i].innerProduct(Rvec[j])=0, starting from low indices, i.e. lowest are least modified
    ///<br> Unique phases: phase=0 at larges component of Rvec[j]
    static std::string orthonormalize(std::vector<Coefficients*> &Dual, std::vector<Coefficients*> &Rvec,
                                      double EliminateSingular=-1. /** set >=0 to eliminate singular vectors, rather than abort */);
    static void normalize(std::vector<Coefficients*> &Dual, std::vector<Coefficients*> &Rvec);

    EigenSolverAbstract& compute(const OperatorAbstract * Op, const OperatorAbstract* Ovr=0); ///< solve the eigenproblem
    bool verify(double Epsilon=1.e-10) const; ///< verify the eigenvectors
    void normalize(); ///< fix phase and normalize to Rvec[i].adjoint() Ovr Rvec[i]=1, adjust dual
    void orthonormalize(); ///< create orthonormal Rvec/Dual sets
    void phaseFix();  ///< fix the phase to =1 for largest in magnitude coefficient
    void sort(std::string Kind); ///< sort eigenvalue and matching vectors by Kind = SmallReal, SmallAbs, LargeReal, LargeAbs etc., see tools::sortKey for all options
    void select(std::string Kind); ///< eigenvalue post-selection (and sorting)

    void parallel(bool Parallel){_serial=not Parallel;} ///< true: use parallelism on block-diagonal operator
    const OperatorAbstract* oper() const {return _op;}
    const OperatorAbstract* ovrl() const {return _ovr;}
};

#endif //EIGEN_SOLVER_ABSTRACT_H
