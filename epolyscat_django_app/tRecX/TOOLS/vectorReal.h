// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef VECTORREAL_H
#define VECTORREAL_H

#include <vector>

/// \ingroup Linalg
/// \brief adds linear space operations to std::vector<double >
class VectorReal:public std::vector<double>
{
public:
    ~VectorReal(){}
    VectorReal(const vector<double> & B):std::vector<double>(B.begin(),B.end()){}
    VectorReal(unsigned int Size=0):std::vector<double>(Size){}
    VectorReal(const VectorReal & B, int Begin, int End):std::vector<double>(End-Begin){for(int k=0;k<End-Begin;k++)data()[k]=B[Begin+k];}
    VectorReal & operator+=(const VectorReal & Other);
    VectorReal & operator-=(const VectorReal & Other);
    VectorReal operator-(const VectorReal & Other) const;

    VectorReal & operator*=(double A);
    VectorReal operator*(double A) const;
    VectorReal & axpy(double A, const VectorReal & X);
    double dot(const VectorReal & X) const;

    double normSqu() const;
    double maxAbsVal() const;

    VectorReal & purge(double Eps=1.e-12); ///< set near-zeros to =0

    static void Test();
};

#endif // VECTORREAL_H
