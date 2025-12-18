// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorFloor.h"
#include "basicDisc.h"
#include "readInput.h"
#include "useMatrix.h"
#include "operatorZD.h"

class OperatorFloorGP : public OperatorFloor
{
    std::vector<double> _basSq;   // norm of basis values at DVR points
    std::vector<std::complex<double>> _GP; // product of basis norms squared with coefficient norms
    std::vector<std::complex<double>> _GP1;
    std::vector<double> _dvrPoints;
    std::vector<double> _dvrWeights;
    std::vector<double> _nodeVal;
    int _nBeg;
    int _bsize;
    int _nSib;
    const Index* _idx;
    void updateNonLin(double time, Coefficients* C);
public:
    OperatorFloorGP(const std::string Name, const std::string Def, const Index *IIndex, const Index *JIndex);
    OperatorFloorGP(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf);
    OperatorFloorGP(std::string Pot, const Index* IIndex, const Index* JIndex, std::complex<double> Multiplier);
    
    void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
              const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const;
    
    void pack(std::vector<int> & Info,std::vector<std::complex<double> > &Buf) const;
    
    static void store(bool store){_store=store;}
    static bool _store;
    static bool _iterations;
};

