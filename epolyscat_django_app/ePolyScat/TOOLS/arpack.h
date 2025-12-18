// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ARPACK_H
#define ARPACK_H

#include <string>
#include <map>
#include <complex>
#include <vector>

/// @brief abstract base class for arpack eigensolver
class Arpack{

    /// @brief basic interface to the fortran77 arpack routines (hidden from outside Arpack)
    class InterfaceF77
    {
        static std::map<std::string,std::string>whichList; ///< which part of the spectrum
        double _tolerance;
        unsigned int _maxIter;
        bool _seriel;

    public:
        InterfaceF77():_tolerance(1.e-12),_maxIter(1000){whichListSet();}
        InterfaceF77( double Tolerance,  /**< accuracy control parameter */
                         int MaxIter    /**< maximal number of arnoldi iterations */,
                         bool Seriel
                         ): _tolerance(Tolerance),_maxIter(MaxIter),_seriel(Seriel){whichListSet();}
        void iter(Arpack *, int K,   /**< desired number of eigenvectors */
                   std::string Which, /**< wich part of the spectrum */
                   bool Restart,      /**< true: iteration starts from first vector in  Rvec */
                   bool computeVectors, bool Parallel=true);


        void whichListSet();
    };
    InterfaceF77 _I_F77;  ///< interface to F77 Arpack

protected:
    unsigned int _lvec;                     ///< number of coefficients a single vector
    std::vector<std::complex<double> > eval;///< eigenvalues
    std::vector<std::complex<double> > rvec;///< storage for right eigenvectors

    /// run the iteration through the reverse-communication InterfaceF77
    void eigenIter(unsigned int Nvec, std::string Which,bool Restart,bool ComputeVectors);

    /// apply operator to vector of length _lvec
    virtual void apply(const std::complex<double> * X, std::complex<double> * Y)=0;

public:
    Arpack(double Tolerance,unsigned int MaxIter,bool Seriel=false)
        :_I_F77(InterfaceF77(Tolerance,MaxIter,Seriel)){}

    /// return eigenvalues as std::vector
    void eigenValues(
            std::vector<std::complex<double> > & Eval,
            unsigned int Nvec=1, ///< number of eigenvectors, =0: all
            const std::string & Which="SmallReal"
            );

    /// return eigenvalues and (right) eigenvectors as std::vector's
    void eigenStdVec(
            std::vector<std::complex<double> > & Eval,
            std::vector<std::vector<std::complex<double> > >& Rvec,
            unsigned int Nvec=1, ///< number of eigenvalues
            const std::string & Which="SmallReal",
            bool Restart=false
            );

    /// verify eigenvectors by direct application and get accuracies
    void verify();
};

/// @brief Example Arpack class using only standard vectors
class ArpackStandard:private Arpack {
    std::vector<std::complex<double> >mat; // this is how this particular instance stores the matrix

    /// only this must be implmented for a given instance
    void apply(const std::complex<double> *X, std::complex<double> * Y);
public:
    ArpackStandard(std::vector<std::complex<double> > Mat,double Tolerance,unsigned int MaxIter)
        :Arpack(Tolerance,MaxIter,true),mat(Mat)
    {_lvec=(unsigned int) sqrt(double(Mat.size())+1);}


    static void test();    ///< basic test
};

#endif // ARPACK_H
