// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARALLELCONTINUITY_H
#define PARALLELCONTINUITY_H

#include <string>
#include <vector>
#include <complex>
#include <map>

#include "index.h"


//class Index;
class Coefficients;

/** \ingroup Parallelization */
/// for imposing continuity
class ParallelContinuity
{

private:
    std::vector<std::vector<std::complex<double>*> > locCur; // pointers to local margin entries
    std::vector<std::vector<std::complex<double>*> > locNei; // pointers to local margin entries

    std::vector<std::vector<std::vector<std::vector<std::complex<double>*> > > > pMarg; // pointers to remote margin entries
    std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > > sendBuf;
    std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > > recvBuf;
    void setIndex(std::string HashC, const Index* I, unsigned int Dimension, int AtBoundary);
    void addNeighbor(std::string HashC, Index * Low, Index * Upp, unsigned int Dimension, unsigned int Level=0);

    // info about lower neighbor of a given index
    class Marg{
        const Index* _cur;  // index for current dimension (at or below Coefficents floor)
        const Index* _nei;  // neighbor of Marg::cur (at or below floor)
        int pCur,pNei;       // positions in floor
        int iCur,iNei;       // positions in index
        int _size; // current size
    public:
        const Index* fNei;  // floor Index hosting Marg::nei (used to retrieve neighbor Coefficients)
        Marg(const Index * Cur,const Index * Nei, const Index * FNei);
        Marg(const Index * PCur,int NCur, const Index * PNei, int NNei, const Index * FNei);
        int curPosFloor();//{return pCur;}
        int neiPosFloor();//{return pNei;}
        int curPosIndex();//{return iCur;}
        int neiPosIndex();//{return iNei;}
        int curSize();//{return _size;}
    };
    // vector of neighbor maps for all indices for directions d<neighbor.size()
    static std::map<std::string, std::vector<std::map<std::string,std::vector<Marg> > > >neighbor;

public:
    ParallelContinuity(){}
    ParallelContinuity(Coefficients *C, int AtBoundary=-1);
    void apply(Coefficients * C, double Scal);

    void setMargin(std::complex<double> Val); ///< set the margin values = Val (non-margin sites remain unused)
    void margin(ParallelContinuity* Marg) const; ///< write present margin into Marg (non-margin sites remain unchanged)
    void halfDiffMargin(ParallelContinuity* Marg) const; ///< write half-difference of margins into Marg  (non-margin sites remain unchanged)
};

#endif // CONTINUITYFLAT_H
