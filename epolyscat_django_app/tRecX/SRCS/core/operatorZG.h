// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORZG_H
#define OPERATORZG_H

#include <string>
#include <vector>
#include <complex>

#include "operatorFloor.h"

class UseMatrix;
//class Operator;

/** \ingroup OperatorFloors */
/// complex general matrix
class OperatorZG: public OperatorFloor
{
    friend class OperatorZGSblock;
    friend class OperatorView;

    static void matmul(const std::complex<double>& A,std::complex<double>*M,std::complex<double> *X,unsigned int SizX, std::complex<double>*Y,unsigned int SizY);
protected:
    void axpy(const std::complex<double> &Alfa, const std::complex<double>* X, unsigned int SizX,
              const std::complex<double> &Beta, std::complex<double>* Y, unsigned int SizY) const;
    void axpyTranspose(const std::complex<double> &Alfa, const std::complex<double>* X, unsigned int SizX,
                       const std::complex<double> &Beta, std::complex<double>* Y, unsigned int SizY) const;
    void axpyTranspose(std::complex<double> Alfa, const std::vector<std::complex<double> > & X,
                       std::complex<double> Beta, std::vector<std::complex<double> > & Y) const{
        axpyTranspose(Alfa,const_cast<std::vector<std::complex<double> >* >(&X)->data(),X.size(),Beta,Y.data(),Y.size());
    }
    void construct(const Eigen::MatrixXcd & Mat, std::string Kind);
public:
    OperatorZG():OperatorFloor(0,0,"ZG"){}
    OperatorZG(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf);
    OperatorZG(const std::vector<std::complex<double> > &Dat,unsigned int Rows, unsigned int Cols, std::string Kind);
    OperatorZG(const UseMatrix *Mat, std::string Kind);
    void pack(std::vector<int> & Info, std::vector<std::complex<double> > &Buf) const;

    void uniqueData();
        
    long applyCount() const;
};

#endif // OPERATORZG_H
