// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DENSITYMATRIX1_H
#define DENSITYMATRIX1_H

#include <memory>
#include <map>
#include <vector>
#include <complex>

#include "qtEigenDense.h"

class BasisOrbitalNumerical;
class BasisIntegrable;
class Index;

///@brief Multi-channel single particle density matrix for an expansion in spherical harmonics
///
/// All parts at a given pair of points (r,s) will be computed as needed:
///
///        rho{ca,cb}[mi,mj,li,lj,i,j](r,s) = Y[li,mi] phi[i](r) rho{cacb}[ij] Y^*[lj,mj] phi^*[j](s)
///
/// where phi[i] is assumed to be in a DVR basis
///
/// designed for accumulating all values at one single pair of radial points (r,s)
class DensityMatrix1{

    Eigen::MatrixXi _cacb; // unique numbering of non-zero blocks of density matrix
    std::vector<Eigen::MatrixXcd> _iSvd,_jSvd; // SVD's for densities all channel pairs ca,cb, labelled by _cacb(ca,cb)

    // current points rho(r,s)
    int _s,_r;

    // lower and upper boundaries of m-expansion
    int _imMin,_imMax;
    int _jmMin,_jmMax;
    // _mlIb[mi][li](i,b)...b'th value of orbital i in range of present IdxB
    std::map<int,std::map<int,Eigen::MatrixXcd> >  _mlIs,_mlJr;
public:
    class M;
    class L;
    class LL{
        friend class DensityMatrix1;
        const M* _densM;
        const Eigen::MatrixXcd *_Is,*_Jr;
        /// density matrix parts at point a,b for given Mi,Mj,Li,Lj
        LL(const L & DensL,int Lj);
    public:
        /// rho[Mi,Mj,Li,Lj,CaCb] at a,b for all channel pairs CaCb
        std::vector<std::complex<double> > densC;
    };

    class L{
        friend class DensityMatrix1;
        int _llMax;
        std::vector<int> _listLj;
        const M* _densM;
        const Eigen::MatrixXcd *_Is;
        const std::map<int,Eigen::MatrixXcd> *_lJr;
        /// density matrix parts at point a,b for given Mi,Mj,Li
        L(M & DensM,int Li);
        std::vector<std::shared_ptr<LL>> _ll;///< rho[Mi,Mj,Li,Lj] at a,b
    public:
        LL* operator() (int Lj); ///< rho[Mi,Mj,Li,Lj] at a,b for all Ca,Cb
        const std::vector<int> & listLj(){return _listLj;}
    };

    class M{
        friend class DensityMatrix1;
        const DensityMatrix1* _dens;
        // lower and upper boundaries of l-expansion
        int _ilMin,_ilMax;
        int _jlMin,_jlMax;
        std::vector<int> _listLi,_listLj;
        const std::map<int,Eigen::MatrixXcd> *_lIs,*_lJr;
        /// density matrix parts at point a,b for given Mi,Mj
        M(DensityMatrix1 &Dens, int Mi, int Mj);
    public:
        L* operator() (int Li); ///< rho[Mi,Mj,Li] at a,b for all Lj,Ca,Cb

        int ilMax() const {return _ilMax;}
        int ilMin() const {return _ilMin;}
        int jlMax() const {return _jlMax;}
        int jlMin() const {return _jlMin;}
        const std::vector<int> & listLi() const {return _listLi;}
        const std::vector<int> & listLj() const {return _listLj;}
        bool isZero() const {return _ilMin>_ilMax or  _jlMin>_jlMax ;}

    private:
        void reset();
        std::map<int,std::shared_ptr<L> > _l;
    };
public:
    DensityMatrix1(const std::vector<std::vector<Eigen::MatrixXcd> > &RhoIJ);
    /// set radial patch defined by BasA and BasB, extract all angular parts from respective OrbI,OrbJ
    void patch(const BasisIntegrable &BasA, const BasisIntegrable &BasB,
               const BasisOrbitalNumerical &OrbI, const BasisOrbitalNumerical &OrbJ);
    M* operator() (int Mi,int Mj);  ///< rho[Mi,Mj] at a,b for all li,lj
    void setPoint(int R, int S){if(_s!=S or _r!=R)reset();_s=S;_r=R;} ///< set new point, discard all stored values

    int nChan() const {return _cacb.cols();}
    int cacb(int ChanA,int ChanB) const {return _cacb(ChanA,ChanB);} ///< consecutive numbering of channel pairs
    int cacbMax() const {return _cacb.lpNorm<Eigen::Infinity>();}

    int imMax() const {return _imMax;} ///< maximal im in current patch
    int imMin() const {return _imMin;} ///< minimal im in current patch
    int jmMax() const {return _jmMax;} ///< maximal im in current patch
    int jmMin() const {return _jmMin;} ///< minimal jm in current patch
private:
    std::map<int,std::map<int,std::shared_ptr<M> >  > _m;
    void reset();
};

#endif // DENSITYMATRIX1_H
